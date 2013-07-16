//
//  Population.cpp
//

#include "config.h"

#include <vector>
#include <algorithm>
#include <iomanip>
#include <RcppArmadillo.h>
#include "SynchronizedUnifGenerator__0__1.h"

#include "Population.h"

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif

// The threads need about
// #define REQUIRED_STACK_SIZE 5242880

using namespace Rcpp;

/*
 * R user interrupt handling helpers
 */
static inline void check_interrupt_impl(void* /*dummy*/) {
	R_CheckUserInterrupt();
}

inline bool check_interrupt() {
	return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}

Population::Population(const Control &ctrl, ::Evaluator &evaluator, SynchronizedUnifGenerator__0__1& unifGen) : ctrl(ctrl), evaluator(evaluator), unifGen(unifGen) {
	// initialize original population (generation 0) totally randomly
	this->currentGeneration.reserve(this->ctrl.populationSize);
	this->nextGeneration.reserve(this->ctrl.populationSize);
	this->currentGenFitnessMap.reserve(this->ctrl.populationSize);
	this->nextGenFitnessMap.reserve(this->ctrl.populationSize);

	this->minCurrentGenFitness = 0.0;
	this->minEliteFitness = 0.0;

#ifdef HAVE_PTHREAD_H
	/*
	 * init mutex for "queue" synchronization
	 * init mutex for Rcout ---- printMutex
	 * init mutex for updating --- updateMutex
	 */
//	if(this->ctrl.numThreads > 1) {
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
//	}
#endif
}

Population::~Population() {
#ifdef HAVE_PTHREAD_H
	this->cleanupPthread();
#endif
}

inline Chromosome& Population::drawChromosomeFromCurrentGeneration(SynchronizedUnifGenerator__0__1& unifGen) {
	int imin = 0, imax = this->ctrl.populationSize - 1;
	int imid = 0;
	
	// Draw a random number between 0 and cumulative sum of all fitness values
	double rand = unifGen() * this->sumCurrentGenFitness;

	// Search for the chromosome whose fitness range surrounds the random number
	while(imin <= imax) {
		imid = (imax + imin) / 2;

		if(this->currentGenFitnessMap[imid] > rand) {
			imax = imid - 1;
		} else {
			imin = imid + 1;
		}
	}

	return this->currentGeneration[imid];
}

inline void Population::addChromosomeToElite(Chromosome& ch) {
	// Add chromosome to elite if better than the worst elite-chromosome
	if(this->ctrl.elitism > 0 && (ch.getFitness() > this->minEliteFitness || this->elite.size() < this->ctrl.elitism)) {
		if(this->elite.size() >= this->ctrl.elitism) {
			this->elite.erase(this->elite.begin());
		}

		this->elite.insert(ch); // The insert copies the chromosome!

		this->minEliteFitness = this->elite.begin()->getFitness();

		if(this->ctrl.verbosity >= MORE_VERBOSE) {
			Rcout << "Adding chromosome to elite. New minimum fitness for elite is " << this->minEliteFitness << std::endl;
		}
	}
}

inline void Population::transformCurrentGenFitnessMap() {
	this->sumCurrentGenFitness = 0.0;
	std::vector<double>::iterator fitnessMapIt;
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity >= DEBUG_VERBOSE) {
		Rcout << "Fitness map: ";
	}
#endif
	
	for(fitnessMapIt = this->currentGenFitnessMap.begin(); fitnessMapIt != this->currentGenFitnessMap.end(); ++fitnessMapIt) {
		this->sumCurrentGenFitness += (*fitnessMapIt) - this->minCurrentGenFitness;
		(*fitnessMapIt) = this->sumCurrentGenFitness;
		
#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.verbosity >= DEBUG_VERBOSE) {
			Rcout << this->sumCurrentGenFitness << " | ";
		}
#endif
	}
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity >= DEBUG_VERBOSE) {
		Rcout << std::endl;
	}
#endif
}

void Population::mate(uint16_t numMatingCouples, ::Evaluator& evaluator, SynchronizedUnifGenerator__0__1& unifGen, bool checkUserInterrupt) {
#ifdef HAVE_PTHREAD_H
	int pthreadRC = 1;
#endif

	for(; numMatingCouples != 0; --numMatingCouples) {
		Chromosome tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(unifGen);
		Chromosome tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(unifGen);
		
		std::vector<Chromosome> children = tmpChromosome1.mateWith(tmpChromosome2, unifGen);
		uint16_t matingTries = 0;
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
		
#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.verbosity >= DEBUG_VERBOSE) {
			Rcout << "Mating chromosomes " << std::endl << tmpChromosome1 << " and" << std::endl << tmpChromosome2 << std::endl
			<< "with minimal fitness " << minParentFitness << std::endl
			<< "First two proposals have fitness " << children[0].getFitness() << " / " << children[1].getFitness() << std::endl;
		}
#endif
		
		// At least the first child should be better than the worse parent
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
			
#ifdef ENABLE_DEBUG_VERBOSITY
			if(this->ctrl.verbosity >= DEBUG_VERBOSE) {
				Rcout << "Proposed children have fitness: " << proposalChildren[0].getFitness() << " / " << proposalChildren[1].getFitness() << std::endl
				<< "Currently selected children have fitness: " << children[0].getFitness() << " / " << children[1].getFitness() << std::endl;
			}
#endif
		}
		
		if(children[0].mutate(unifGen) == true) {
			evaluator.evaluate(children[0]);
		}
		if(children[1].mutate(unifGen) == true) {
			evaluator.evaluate(children[1]);
		}
		
		if(this->ctrl.verbosity >= MORE_VERBOSE) {
#ifdef HAVE_PTHREAD_H
			pthreadRC = pthread_mutex_lock(&this->printMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
#endif
			this->printChromosomeFitness(Rcout, children[0]);
			this->printChromosomeFitness(Rcout, children[1]);
#ifdef HAVE_PTHREAD_H
			pthreadRC = pthread_mutex_unlock(&this->printMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
#endif
		}

		/*
		 * Update fitness map and add children to the generation
		 */		
#ifdef HAVE_PTHREAD_H
		pthreadRC = pthread_mutex_lock(&this->updateMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
#endif
		this->nextGenFitnessMap.push_back(children[0].getFitness());
		this->nextGenFitnessMap.push_back(children[1].getFitness());
		
		if(children[0].getFitness() < this->minCurrentGenFitness) {
			this->minCurrentGenFitness = children[0].getFitness();
		}
		
		if(children[1].getFitness() < this->minCurrentGenFitness) {
			this->minCurrentGenFitness = children[1].getFitness();
		}
		
		this->addChromosomeToElite(children[0]);
		this->addChromosomeToElite(children[1]);
		
		this->nextGeneration.push_back(children[0]);
		this->nextGeneration.push_back(children[1]);
#ifdef HAVE_PTHREAD_H
		pthreadRC = pthread_mutex_unlock(&this->updateMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
#endif

		if(checkUserInterrupt == true) {
			if(check_interrupt()) {
				throw InterruptException();
			}
		}
	}
}

void Population::run() {
	int i = 0;
	VariablePositionPopulation varPosPop(this->ctrl.chromosomeSize);
	uint16_t mainThreadMatingCouples = this->ctrl.populationSize / 2;

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
		
		if(check_interrupt()) {
			throw InterruptException();
		}
	}

#ifdef HAVE_PTHREAD_H
	Population::ThreadArgsWrapper* threadArgs;
	uint16_t maxThreadsToSpawn = this->ctrl.numThreads - 1;
	
	int pthreadRC = 1; // return value for of pthread_ functions
	
	if(maxThreadsToSpawn > 0) {
		uint16_t otherThreadsMatingCouples = mainThreadMatingCouples / this->ctrl.numThreads;
		uint16_t remainingMatingCouples = (mainThreadMatingCouples % this->ctrl.numThreads);
		mainThreadMatingCouples = otherThreadsMatingCouples;
				
		threadArgs = new Population::ThreadArgsWrapper[maxThreadsToSpawn];
		this->threads = new pthread_t[maxThreadsToSpawn];
		
		pthread_attr_t threadAttr;

		pthreadRC = pthread_attr_init(&threadAttr);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC);
		if(pthreadRC != 0) {
			throw ThreadingError("Thread attributes could not be initialized");
		}
		
		pthreadRC = pthread_attr_setdetachstate(&threadAttr, PTHREAD_CREATE_JOINABLE);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC);
				
		pthreadRC = pthread_attr_setstacksize(&threadAttr, 1024 * 1024 * 1024 * 5);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC);
		
		if(pthreadRC != 0) {
			throw ThreadingError("Thread attributes could not be modified to make the thread joinable");
		}
		
		for(i = maxThreadsToSpawn - 1; i >= 0; --i) {
			threadArgs[i].numMatingCouples = otherThreadsMatingCouples;
			
			if(remainingMatingCouples-- > 0) {
				++threadArgs[i].numMatingCouples;
			}
			
			threadArgs[i].popObj= this;
			threadArgs[i].unifGen = new SynchronizedUnifGenerator__0__1(UNIF_GENERATOR_BUFFER_SIZE_THREAD);
			threadArgs[i].evalObj = this->evaluator.clone();
			pthreadRC = pthread_create((this->threads + i), &threadAttr, &Population::matingThreadStart, (void *) (threadArgs + i));
			
			if(pthreadRC == 0) {
				++this->actuallySpawnedThreads;
			}
#ifdef ENABLE_DEBUG_VERBOSITY
			else {
				if(this->ctrl.verbosity >= DEBUG_VERBOSE) {
					Rcerr << "Warning: Thread " << i << " could not be created" << std::endl;
				}
			}
#endif
		}
		
		if(this->actuallySpawnedThreads < maxThreadsToSpawn) {
			Rcerr << "Warning: Only " << this->actuallySpawnedThreads << " threads could be spawned" << std::endl;
			mainThreadMatingCouples += otherThreadsMatingCouples * (maxThreadsToSpawn - this->actuallySpawnedThreads);
		}
		
		pthreadRC = pthread_attr_destroy(&threadAttr);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC);

		if(this->ctrl.verbosity >= OFF) {
			Rcout << "Spawned " << this->actuallySpawnedThreads << " threads" << std::endl;
		}
	}
#endif

	bool interrupted = false;

	for(i = this->ctrl.numGenerations; i > 0 && !interrupted; --i) {
		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}

		/*
		 * Transform the fitness map of the current generation to start at 0
		 * and have cumulative values
		 */
		this->transformCurrentGenFitnessMap();
		this->minCurrentGenFitness = 0.0;

#ifdef HAVE_PTHREAD_H
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
#endif
		/*
		 * Mate two chromosomes to generate two children that are eventually mutated
		 * To get the same population size, a total of popSize / 2 mating pairs have
		 * to generate 2 children
		 *
		 */
		try {
			this->mate(mainThreadMatingCouples, this->evaluator, this->unifGen, true);
		} catch (InterruptException) {
			interrupted = true;
		}

#ifdef HAVE_PTHREAD_H		
		this->finishedMating();
		this->waitForAllThreadsToFinishMating();
#endif
		// Housekeeping
		this->currentGeneration = this->nextGeneration;
		this->nextGeneration.clear();
		
		this->currentGenFitnessMap = this->nextGenFitnessMap;
		this->nextGenFitnessMap.clear();
	}

#ifdef HAVE_PTHREAD_H
	/*
	 * signal threads to end
	 */
	if(maxThreadsToSpawn > 0) {
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
	}
#endif
	
	if(interrupted) {
		throw InterruptException();
	}
}

SortedChromosomes Population::getResult() const {
	SortedChromosomes result(this->currentGeneration.begin(), this->currentGeneration.end());

	if(this->ctrl.elitism > 0 && this->elite.size() > 0) {
		result.insert(this->elite.begin(), this->elite.end());
	}

	return result;
}

inline std::ostream& Population::printChromosomeFitness(std::ostream &os, Chromosome &ch) {
	os << ch << TAB_DELIMITER << std::fixed
		<< std::setw(WIDTH) << std::setprecision(PRECISION) << ch.getFitness() << std::endl;

	return os;
}


#ifdef HAVE_PTHREAD_H
void* Population::matingThreadStart(void* obj) {
	ThreadArgsWrapper* args = static_cast<ThreadArgsWrapper*>(obj);
	args->evalObj->setUnifGenerator(args->unifGen);
	args->popObj->runMating(args->numMatingCouples, *args->evalObj, *(args->unifGen));
	return NULL;
}

inline void Population::cleanupPthread() {
//	if(this->ctrl.numThreads > 1) {
		pthread_mutex_destroy(&this->printMutex);
		pthread_mutex_destroy(&this->updateMutex);
		pthread_mutex_destroy(&this->syncMutex);
		pthread_cond_destroy(&this->startMatingCond);
		pthread_cond_destroy(&this->allThreadsFinishedMatingCond);
//	}
}

inline void Population::cancelAllThreads() {
	if(this->actuallySpawnedThreads > 0) {
		/*
		 * First try to softly shut down the threads
		 */
		int pthreadRC = pthread_mutex_lock(&this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		this->startMating = true;
		this->killThreads = true;
		
		pthreadRC = pthread_cond_broadcast(&this->startMatingCond);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		pthreadRC = pthread_mutex_unlock(&this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		/*
		 * Now forcefully cancel the threads
		 */
		for(int i = this->ctrl.numThreads - 2; i >= 0; --i) {
			/*
			 * If the thread was never created (i.e. pthread_create failed) the call will return an
			 * error code, but it will not block the thread!
			 */
			pthreadRC = pthread_cancel(this->threads[i]);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		}
	}
	this->cleanupPthread();
}

void Population::runMating(uint16_t numMatingCouples, ::Evaluator& evaluator, SynchronizedUnifGenerator__0__1& unifGen) {
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
		this->mate(numMatingCouples, evaluator, unifGen, false);
		
		/*
		 * Signal that the thread has finished mating
		 */
		this->finishedMating();
		this->waitForAllThreadsToFinishMating();
	}
}

inline void Population::waitForAllThreadsToFinishMating() {
	int pthreadRC = pthread_mutex_lock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	while(this->allThreadsFinishedMating == false) {
		pthreadRC = pthread_cond_wait(&this->allThreadsFinishedMatingCond, &this->syncMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	}
	
	pthreadRC = pthread_mutex_unlock(&this->syncMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
}

inline void Population::finishedMating() {
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


