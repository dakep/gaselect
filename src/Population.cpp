//
//  Population.cpp
//

#include "config.h"

#include <vector>
#include <algorithm>
#include <iomanip>
#include <RcppArmadillo.h>
#include "UnifGenerator__0__1.h"

#include "Population.h"

#ifdef HAVE_PTHREADS
#include <pthread.h>
#endif

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

Population::Population(const Control &ctrl, const ::Evaluator &evaluator) : unifGen(), ctrl(ctrl), evaluator(&evaluator) {
	// initialize original population (generation 0) totally randomly
	this->currentGeneration.reserve(this->ctrl.populationSize);
	this->nextGeneration.reserve(this->ctrl.populationSize);
	this->currentGenFitnessMap.reserve(this->ctrl.populationSize);
	this->nextGenFitnessMap.reserve(this->ctrl.populationSize);

	this->minCurrentGenFitness = 0.0;
	this->minEliteFitness = 0.0;
}

Population::~Population() {
//	SortedChromosomes::iterator rmEliteIter = this->elite.begin();
//
//	// Delete all elite chromosomes from stack
//	for(; rmEliteIter != this->elite.end(); ++rmEliteIter) {
//		delete (*rmEliteIter);
//	}
//
//	this->cleanCurrentGeneration();
}

inline Chromosome& Population::drawChromosomeFromCurrentGeneration() {
	int imin = 0, imax = this->ctrl.populationSize - 1;
	int imid = 0;
	
	// Draw a random number between 0 and cumulative sum of all fitness values
	double rand = this->unifGen() * this->sumCurrentGenFitness;

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
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << "Fitness map: ";
	}
#endif
	
	for(fitnessMapIt = this->currentGenFitnessMap.begin(); fitnessMapIt != this->currentGenFitnessMap.end(); ++fitnessMapIt) {
		this->sumCurrentGenFitness += (*fitnessMapIt) - this->minCurrentGenFitness;
		(*fitnessMapIt) = this->sumCurrentGenFitness;
		
#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.verbosity == DEBUG_VERBOSE) {
			Rcout << this->sumCurrentGenFitness << " | ";
		}
#endif
	}
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << std::endl;
	}
#endif
}


inline void Population::mate(uint16_t numMatingCouples) {
	for(; numMatingCouples != 0; --numMatingCouples) {
		Chromosome tmpChromosome1 = this->drawChromosomeFromCurrentGeneration();
		Chromosome tmpChromosome2 = this->drawChromosomeFromCurrentGeneration();
		
		std::vector<Chromosome> children = tmpChromosome1.mateWith(tmpChromosome2);
		uint16_t matingTries = 0;
		double minParentFitness = ((tmpChromosome1.getFitness() > tmpChromosome2.getFitness()) ? tmpChromosome1.getFitness() : tmpChromosome2.getFitness());
		
		this->evaluator->evaluate(children[0]);
		this->evaluator->evaluate(children[1]);
		// Make sure the first child is "better" than the second child
		if(children[0].getFitness() < children[1].getFitness()) {
			std::swap(children[1], children[0]);
		}
		
		
#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.verbosity == DEBUG_VERBOSE) {
			Rcout << "Mating chromosomes " << std::endl << tmpChromosome1 << " and" << std::endl << tmpChromosome2 << std::endl
			<< "with minimal fitness " << minParentFitness << std::endl
			<< "First two proposals have fitness " << children[0].getFitness() << " / " << children[1].getFitness() << std::endl;
		}
#endif
		
		// At least the first child should be better than the worse parent
		while((children[0].getFitness() < minParentFitness) && (++matingTries < this->ctrl.maxMatingTries)) {
			std::vector<Chromosome> proposalChildren = tmpChromosome1.mateWith(tmpChromosome2);
			
			/*
			 * After mating a chromosome may have no variables at all, so we need to check if the variable count is
			 * greater than 0, otherwise the evaluation step would fail
			 */			
			if((proposalChildren[0].getVariableCount() > 0) && (this->evaluator->evaluate(proposalChildren[0]) > children[1].getFitness())) { // better as 2nd child
				if(proposalChildren[0].getFitness() > children[0].getFitness()) { // even better as 1st child
					children[1] = children[0];
					children[0] = proposalChildren[0];
				} else {
					children[1] = proposalChildren[0];
				}
			}
			
			// Check 2nd new child
			if((proposalChildren[1].getVariableCount() > 0) && (this->evaluator->evaluate(proposalChildren[1]) > children[1].getFitness())) { // better as 2nd child
				if(proposalChildren[1].getFitness() > children[0].getFitness()) { // even better as 1st child
					children[1] = children[0];
					children[0] = proposalChildren[1];
				} else {
					children[1] = proposalChildren[1];
				}
			}
			
#ifdef ENABLE_DEBUG_VERBOSITY
			if(this->ctrl.verbosity == DEBUG_VERBOSE) {
				Rcout << "Proposed children have fitness: " << proposalChildren[0].getFitness() << " / " << proposalChildren[1].getFitness() << std::endl
				<< "Currently selected children have fitness: " << children[0].getFitness() << " / " << children[1].getFitness() << std::endl;
			}
#endif
		}
		
		if(children[0].mutate() == true) {
			this->evaluator->evaluate(children[0]);
		}
		if(children[1].mutate() == true) {
			this->evaluator->evaluate(children[1]);
		}
		
		this->nextGenFitnessMap.push_back(children[0].getFitness());
		this->nextGenFitnessMap.push_back(children[1].getFitness());
		
		if(children[0].getFitness() < this->minCurrentGenFitness) {
			this->minCurrentGenFitness = children[0].getFitness();
		}
		
		if(children[1].getFitness() < this->minCurrentGenFitness) {
			this->minCurrentGenFitness = children[1].getFitness();
		}
		
		if(this->ctrl.verbosity >= MORE_VERBOSE) {
			this->printChromosomeFitness(Rcout, children[0]);
			this->printChromosomeFitness(Rcout, children[1]);
		}
		
		this->addChromosomeToElite(children[0]);
		this->addChromosomeToElite(children[1]);
		
		this->nextGeneration.push_back(children[0]);
		this->nextGeneration.push_back(children[1]);
		
		if(check_interrupt()) {
			throw InterruptException();
		}
	}
}

void Population::run() {
	int i = 0;
	int popSizeHalf = this->ctrl.populationSize / 2;
	VariablePositionPopulation varPosPop(this->ctrl.chromosomeSize);

	if(this->ctrl.verbosity > OFF) {
		Rcout << "Generating initial population" << std::endl;
	}

	for(i = this->ctrl.populationSize; i > 0; --i) {
		Chromosome tmpChromosome(this->ctrl, varPosPop);

		this->currentGenFitnessMap.push_back(this->evaluator->evaluate(tmpChromosome));
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

#ifdef HAVE_PTHREADS
	/*
	 * init mutex for "queue" synch
	 * init mutex for printing
	 * create numThreads threads pthread_create(&tid[i], NULL, Population::runMating, this)
	 *
	 */
#endif
	
	for(i = this->ctrl.numGenerations; i > 0; --i) {
		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}

		/*
		 * Transform the fitness map of the current generation to start at 0
		 * and have cumulative values
		 */
		this->transformCurrentGenFitnessMap();
		this->minCurrentGenFitness = 0.0;
		
		/*
		 * Mate two chromosomes to generate two children that are eventually mutated
		 * To get the same population size, a total of popSize / 2 mating pairs have
		 * to generate 2 children
		 *
		 */
		this->mate(popSizeHalf);

		// Housekeeping
		this->currentGeneration = this->nextGeneration;
		this->nextGeneration.clear();

		this->currentGenFitnessMap = this->nextGenFitnessMap;
		this->nextGenFitnessMap.clear();
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
