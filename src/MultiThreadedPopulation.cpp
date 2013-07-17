//
//  MultiThreadedPopulation.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//
#include "config.h"

#define HAVE_PTHREAD_H 1

#ifdef HAVE_PTHREAD_H

#include <exception>
#include <vector>
#include <algorithm>
#include <RcppArmadillo.h>
#include <pthread.h>
#include <errno.h>

#include "UnifGenerator_0_1.h"
#include "MultiThreadedPopulation.h"

using namespace Rcpp;

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->ctrl.verbosity >= DEBUG_VERBOSE) { expr; }
#else
#define IF_DEBUG(expr)
#endif

/*
 * R user interrupt handling helpers
 */
static inline void check_interrupt_impl(void* /*dummy*/) {
	R_CheckUserInterrupt();
}

inline bool check_interrupt() {
	return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

MultiThreadedPopulation::MultiThreadedPopulation(const Control &ctrl, ::Evaluator &evaluator, SynchronizedUnifGenerator_0_1& unifGen) : Population(ctrl, evaluator), unifGen(unifGen), nextGenFitnessMap(ctrl.populationSize, 0.0) {
	// initialize original population (generation 0) totally randomly
	if(this->ctrl.numThreads <= 1) {
		throw new std::logic_error("This population should only be used if multiple threads are requested");
	}
	
	this->nextGeneration.reserve(this->ctrl.populationSize);
	
	this->minCurrentGenFitness = 0.0;
	
	int pthreadRC = pthread_mutex_init(&this->printMutex, NULL);
	if(pthreadRC != 0) {
		throw ThreadingError("Mutex for printing could not be initialized");
	}
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	
	pthreadRC = pthread_mutex_init(&this->updateMutex, NULL);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	if(pthreadRC != 0) {
		throw ThreadingError("Mutex for updating could not be initialized");
	}
	
	pthreadRC = pthread_mutex_init(&this->syncMutex, NULL);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	if(pthreadRC != 0) {
		throw ThreadingError("Mutex for synchronization could not be initialized");
	}
	
	pthreadRC = pthread_cond_init(&this->startMatingCond, NULL);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	if(pthreadRC != 0) {
		throw ThreadingError("Condition for synchronization (start mating) could not be initialized");
	}
	
	pthreadRC = pthread_cond_init(&this->allThreadsFinishedMatingCond, NULL);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	if(pthreadRC != 0) {
		throw ThreadingError("Condition for synchronization (finished mating) could not be initialized");
	}
	
	this->startMating = false;
	this->allThreadsFinishedMating = false;
	this->killThreads = false;
	
	this->actuallySpawnedThreads = 0;
	this->numThreadsFinishedMating = 0;
}

MultiThreadedPopulation::~MultiThreadedPopulation() {
	pthread_mutex_destroy(&this->printMutex);
	pthread_mutex_destroy(&this->updateMutex);
	pthread_mutex_destroy(&this->syncMutex);
	pthread_cond_destroy(&this->startMatingCond);
	pthread_cond_destroy(&this->allThreadsFinishedMatingCond);
}

inline void MultiThreadedPopulation::transformCurrentGenFitnessMap() {
	this->sumCurrentGenFitness = 0.0;
	std::vector<double>::iterator fitnessMapIt;
	IF_DEBUG(Rcout << "Fitness map: ")
	
	for(fitnessMapIt = this->currentGenFitnessMap.begin(); fitnessMapIt != this->currentGenFitnessMap.end(); ++fitnessMapIt) {
		this->sumCurrentGenFitness += (*fitnessMapIt) - this->minCurrentGenFitness;
		(*fitnessMapIt) = this->sumCurrentGenFitness;

		IF_DEBUG(Rcout << this->sumCurrentGenFitness << " | ")
	}
	IF_DEBUG(Rcout << std::endl)
}

void MultiThreadedPopulation::mate(uint16_t numMatingCouples, ::Evaluator& evaluator, SynchronizedUnifGenerator_0_1& unifGen, uint16_t offset, bool checkUserInterrupt) {
	int pthreadRC = 1;
	int i = 0;
	
	uint16_t matingTries = 0;
	
	for(; i < numMatingCouples; ++i) {
		Chromosome tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(unifGen() * this->sumCurrentGenFitness);
		Chromosome tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(unifGen() * this->sumCurrentGenFitness);
		
		std::vector<Chromosome> children = tmpChromosome1.mateWith(tmpChromosome2, unifGen);
		double minParentFitness = ((tmpChromosome1.getFitness() > tmpChromosome2.getFitness()) ? tmpChromosome1.getFitness() : tmpChromosome2.getFitness());
		
		/*
		 * If both children have no variables, mate again
		 */
		while(children[0].getVariableCount() == 0 && children[1].getVariableCount() == 0) {
			tmpChromosome1.mateWith(tmpChromosome2, unifGen);
		}
		
		if(children[0].getVariableCount() == 0) {
			children[0] = children[1];
		} else if(children[1].getVariableCount() == 0) {
			children[1] = children[0];
		}
		
		evaluator.evaluate(children[0]);
		evaluator.evaluate(children[1]);
		// Make sure the first child is "better" than the second child
		if(children[0].getFitness() < children[1].getFitness()) {
			std::swap(children[1], children[0]);
		}
		
		IF_DEBUG(
			Rcout << "Mating chromosomes " << std::endl << tmpChromosome1 << " and" << std::endl << tmpChromosome2 << std::endl
				 << "with minimal fitness " << minParentFitness << std::endl
			<< "First two proposals have fitness " << children[0].getFitness() << " / " << children[1].getFitness() << std::endl;
		)
		
		// At least the first child should be better than the worse parent
		matingTries = 0;
		while((children[0].getFitness() < minParentFitness) && (++matingTries < this->ctrl.maxMatingTries)) {
			std::vector<Chromosome> proposalChildren = tmpChromosome1.mateWith(tmpChromosome2, unifGen);
			
			/*
			 * After mating a chromosome may have no variables at all, so we need to check if the variable count is
			 * greater than 0, otherwise the evaluation step would fail
			 */
			if(proposalChildren[0].getVariableCount() > 0) {
				if(evaluator.evaluate(proposalChildren[0]) > children[1].getFitness()) { // better as 2nd child
					if(proposalChildren[0].getFitness() > children[0].getFitness()) { // even better as 1st child
						children[1] = children[0];
						children[0] = proposalChildren[0];
					} else {
						children[1] = proposalChildren[0];
					}
				}
			}
			
			// Check 2nd new child
			if(proposalChildren[1].getVariableCount() > 0) {
				if(evaluator.evaluate(proposalChildren[1]) > children[1].getFitness()) { // better as 2nd child
					if(proposalChildren[1].getFitness() > children[0].getFitness()) { // even better as 1st child
						children[1] = children[0];
						children[0] = proposalChildren[1];
					} else {
						children[1] = proposalChildren[1];
					}
				}
			}
			
			IF_DEBUG(
				Rcout << "Proposed children have fitness: " << proposalChildren[0].getFitness() << " / " << proposalChildren[1].getFitness() << std::endl
					 << "Currently selected children have fitness: " << children[0].getFitness() << " / " << children[1].getFitness() << std::endl;
			)
		}
		
		if(children[0].mutate(unifGen) == true) {
			evaluator.evaluate(children[0]);
		}
		if(children[1].mutate(unifGen) == true) {
			evaluator.evaluate(children[1]);
		}
		
		if(this->ctrl.verbosity >= MORE_VERBOSE) {
			pthreadRC = pthread_mutex_lock(&this->printMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
			
			this->printChromosomeFitness(Rcout, children[0]);
			this->printChromosomeFitness(Rcout, children[1]);
			
			pthreadRC = pthread_mutex_unlock(&this->printMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		}
		
		/*
		 * Update fitness map and add children to the generation
		 *
		 * It is guaranteed that the elements
		 * offset, offset + 1, ..., offset + 2 * i, offet + 2 * i + 1, ..., offset + 2 * numMatingCouples - 1
		 * of the fitness map and the generation are only accessed
		 * by this particular thread, so no locking is neccessary
		 */
		this->nextGenFitnessMap[offset + (2 * i)] = children[0].getFitness();
		this->nextGenFitnessMap[offset + (2 * i) + 1] = children[1].getFitness();

		this->nextGeneration[offset + (2 * i)] = children[0];
		this->nextGeneration[offset + (2 * i) + 1] = children[1];

		/*
		 * The main thread has to check for a user interrupt
		 */
		if(checkUserInterrupt == true) {
			if(check_interrupt()) {
				throw InterruptException();
			}
		}
	}
	
	/*
	 * Acquire lock to recalculate the minimum fitness and add
	 * elements to the elite if appropriate
	 */
	pthreadRC = pthread_mutex_lock(&this->updateMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)

	for(i = offset + (2 * numMatingCouples) - 1; i >= offset; --i) {

		if(this->nextGenFitnessMap[i] < this->minCurrentGenFitness) {
			this->minCurrentGenFitness = this->nextGenFitnessMap[i];
		}
		this->addChromosomeToElite(this->nextGeneration[i]);
	}

	pthreadRC = pthread_mutex_unlock(&this->updateMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
}

void MultiThreadedPopulation::run() {
	int i = 0;
	VariablePositionPopulation varPosPop(this->ctrl.chromosomeSize);
	MultiThreadedPopulation::ThreadArgsWrapper* threadArgs;
	uint16_t maxThreadsToSpawn = this->ctrl.numThreads - 1;
	uint16_t mainThreadMatingCouples = this->ctrl.populationSize / 2;
	uint16_t otherThreadsMatingCouples = mainThreadMatingCouples / this->ctrl.numThreads;
	uint16_t remainingMatingCouples = (mainThreadMatingCouples % this->ctrl.numThreads);
	uint16_t offset = 0;
	pthread_attr_t threadAttr;

	int pthreadRC = 1; // return value for of pthread_ functions

	mainThreadMatingCouples = otherThreadsMatingCouples;
	
	if(this->ctrl.verbosity > OFF) {
		Rcout << "Generating initial population" << std::endl;
	}
	
	for(i = this->ctrl.populationSize; i > 0; --i) {
		Chromosome tmpChromosome(this->ctrl, varPosPop, this->unifGen);
		
		this->currentGenFitnessMap.push_back(this->evaluator.evaluate(tmpChromosome));
		if(tmpChromosome.getFitness() < this->minCurrentGenFitness) {
			this->minCurrentGenFitness = tmpChromosome.getFitness();
		}
		
		if(this->ctrl.verbosity >= MORE_VERBOSE) {
			this->printChromosomeFitness(Rcout, tmpChromosome);
		}
		this->addChromosomeToElite(tmpChromosome);
		this->currentGeneration.push_back(tmpChromosome);
		
		// Full next generation with dummies
		this->nextGeneration.push_back(Chromosome(tmpChromosome, false));
		
		if(check_interrupt()) {
			throw InterruptException();
		}
	}
	
	threadArgs = new MultiThreadedPopulation::ThreadArgsWrapper[maxThreadsToSpawn];
	this->threads = new pthread_t[maxThreadsToSpawn];
	
	
	pthreadRC = pthread_attr_init(&threadAttr);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	if(pthreadRC != 0) {
		throw ThreadingError("Thread attributes could not be initialized");
	}
	
	pthreadRC = pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	
//	pthreadRC = pthread_attr_setstacksize(&threadAttr, 1024 * 1024 * 5);
//	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	
	if(pthreadRC != 0) {
		throw ThreadingError("Thread attributes could not be modified to make the thread joinable");
	}

	for(i = maxThreadsToSpawn - 1; i >= 0; --i) {
		threadArgs[i].numMatingCouples = otherThreadsMatingCouples;
		
		if(remainingMatingCouples > 0) {
			--remainingMatingCouples;
			++threadArgs[i].numMatingCouples;
		}
		threadArgs[i].offset = offset;
		threadArgs[i].popObj = this;
		threadArgs[i].unifGen = new SynchronizedUnifGenerator_0_1(UNIF_GENERATOR_BUFFER_SIZE_THREAD);
		threadArgs[i].evalObj = this->evaluator.clone();

		pthreadRC = pthread_create((this->threads + i), &threadAttr, &MultiThreadedPopulation::matingThreadStart, (void *) (threadArgs + i));
		
		if(pthreadRC == 0) {
			++this->actuallySpawnedThreads;
			offset += 2 * threadArgs[i].numMatingCouples;
		} else {
			mainThreadMatingCouples += threadArgs[i].numMatingCouples;			
			IF_DEBUG(
					 Rcerr << "Warning: Thread " << i << " could not be created: " << strerror(pthreadRC) << std::endl;
			)
		}
	}
		
	if(this->actuallySpawnedThreads < maxThreadsToSpawn) {
		Rcerr << "Warning: Only " << this->actuallySpawnedThreads << " threads could be spawned" << std::endl;
	}
	
	pthreadRC = pthread_attr_destroy(&threadAttr);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC);
	
	if(this->ctrl.verbosity >= OFF) {
		Rcout << "Spawned " << this->actuallySpawnedThreads << " threads" << std::endl;
	}
	
	bool interrupted = false;
	
	for(i = this->ctrl.numGenerations; i > 0 && !interrupted; --i) {
		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}
		
		/*
		 * Transform the fitness map of the current generation to start at 0
		 * and have cumulative values
		 */
		this->transformCurrentGenFitnessMap(); // It is guaranteed that no other thread is running, so it is safe to not use locks!
		this->minCurrentGenFitness = 0.0;
		
		/*
		 * broadcast to all threads to start mating
		 */
		pthreadRC = pthread_mutex_lock(&this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		this->startMating = true;
		
		pthreadRC = pthread_cond_broadcast(&this->startMatingCond);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		pthreadRC = pthread_mutex_unlock(&this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		/*
		 * Mate two chromosomes to generate two children that are eventually mutated
		 * To get the same population size, a total of popSize / 2 mating pairs have
		 * to generate 2 children
		 *
		 */
		try {
			this->mate(mainThreadMatingCouples, this->evaluator, this->unifGen, offset, true);
		} catch (InterruptException) {
			interrupted = true;
		}
		
		this->finishedMating();
		this->waitForAllThreadsToFinishMating();
		
		// Housekeeping
		this->currentGeneration = this->nextGeneration;
		this->currentGenFitnessMap = this->nextGenFitnessMap;
	}
	
	/*
	 * signal threads to end
	 */
	pthreadRC = pthread_mutex_lock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	this->startMating = true;
	this->killThreads = true;
	
	pthreadRC = pthread_cond_broadcast(&this->startMatingCond);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	pthreadRC = pthread_mutex_unlock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	
	for(i = maxThreadsToSpawn - 1; i >= 0; --i) {
		/*
		 * If the thread was never created (i.e. pthread_create failed) the call will return an
		 * error code, but it will not block the thread!
		 */
		pthreadRC = pthread_join(this->threads[i], NULL);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		delete threadArgs[i].unifGen;
		delete threadArgs[i].evalObj;
	}
	
	delete threadArgs;
	
	if(interrupted) {
		throw InterruptException();
	}
}

void* MultiThreadedPopulation::matingThreadStart(void* obj) {
	ThreadArgsWrapper* args = static_cast<ThreadArgsWrapper*>(obj);
	args->evalObj->setUnifGenerator(args->unifGen);
	args->popObj->runMating(args->numMatingCouples, *args->evalObj, *(args->unifGen), args->offset);
	return NULL;
}

void MultiThreadedPopulation::runMating(uint16_t numMatingCouples, ::Evaluator& evaluator, SynchronizedUnifGenerator_0_1& unifGen, uint16_t offset) {
	int pthreadRC = 1;
	while(true) {
		/*
		 * Wait until the thread is started
		 */
		pthreadRC = pthread_mutex_lock(&this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		while(this->startMating == false) {
			pthreadRC = pthread_cond_wait(&this->startMatingCond, &this->syncMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		}
		
		/*
		 * Check if the thread is killed
		 */
		if(this->killThreads == true) {
			pthreadRC = pthread_mutex_unlock(&this->syncMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
			break;
		}
		
		pthreadRC = pthread_mutex_unlock(&this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		/*
		 * Do actual mating
		 */
		this->mate(numMatingCouples, evaluator, unifGen, offset, false);
		
		/*
		 * Signal that the thread has finished mating
		 */
		this->finishedMating();
		this->waitForAllThreadsToFinishMating();
	}
}

inline void MultiThreadedPopulation::waitForAllThreadsToFinishMating() {
	int pthreadRC = pthread_mutex_lock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	while(this->allThreadsFinishedMating == false) {
		pthreadRC = pthread_cond_wait(&this->allThreadsFinishedMatingCond, &this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	}
	
	pthreadRC = pthread_mutex_unlock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
}

inline void MultiThreadedPopulation::finishedMating() {
	int pthreadRC = pthread_mutex_lock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	if(++this->numThreadsFinishedMating > this->actuallySpawnedThreads) { // > because the main thread must finish mating as well
		this->allThreadsFinishedMating = true;
		this->numThreadsFinishedMating = 0;
		this->startMating = false;
		
		pthreadRC = pthread_cond_broadcast(&this->allThreadsFinishedMatingCond);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	} else {
		this->allThreadsFinishedMating = false;
	}
	
	pthreadRC = pthread_mutex_unlock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
}

#endif