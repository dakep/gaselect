//
//  SingleThreadPopulation.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//
#include "config.h"

#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include <RcppArmadillo.h>

#include "RNG.h"
#include "SingleThreadPopulation.h"
#include "ShuffledSet.h"

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

SingleThreadPopulation::SingleThreadPopulation(const Control &ctrl, ::Evaluator &evaluator, const std::vector<uint32_t> &seed) : Population(ctrl, evaluator, seed) {}


void SingleThreadPopulation::run() {
	int i = 0;
	uint16_t j = 0;
	uint8_t matingTries = 0;
	ShuffledSet shuffledSet(this->ctrl.chromosomeSize);
	RNG rng(this->seed);
	
	double sumFitness = 0.0;
	double minFitness = 0.0;
	double minParentFitness = 0.0;
	
	std::vector<Chromosome*> newGeneration;
	
	Chromosome* tmpChromosome1;
	Chromosome* tmpChromosome2;
	std::vector<Chromosome*>::iterator child1It;
	std::vector<Chromosome*>::reverse_iterator child2It;
	Chromosome* proposalChild1 = new Chromosome(this->ctrl, shuffledSet, rng, false);
	Chromosome* proposalChild2 = new Chromosome(this->ctrl, shuffledSet, rng, false);
	uint8_t child1Tries = 0;
	uint8_t child2Tries = 0;
	bool child1Mutated, child2Mutated;
	bool child1Dup, child2Dup;
	std::vector<Chromosome*>::iterator findIt;
	std::vector<Chromosome*>::reverse_iterator revFindIt;
	
	newGeneration.reserve(this->ctrl.populationSize);
	
	if(this->ctrl.verbosity > OFF) {
		Rcout << "Generating initial population" << std::endl;
	}
	
	while(newGeneration.size() < this->ctrl.populationSize) {
		tmpChromosome1 = new Chromosome(this->ctrl, shuffledSet, rng);
		
		/* Check if chromosome is already in the initial population */
		if(std::find_if(newGeneration.begin(), newGeneration.end(), CompChromsomePtr(tmpChromosome1)) == newGeneration.end()) {
			this->currentGenFitnessMap.push_back(this->evaluator.evaluate(*tmpChromosome1));
			if(tmpChromosome1->getFitness() < minFitness) {
				minFitness = tmpChromosome1->getFitness();
			}
			
			if(this->ctrl.verbosity >= MORE_VERBOSE) {
				this->printChromosomeFitness(Rcout, *tmpChromosome1);
			}
			
			this->addChromosomeToElite(*tmpChromosome1);
			
			newGeneration.push_back(tmpChromosome1);
			this->currentGeneration.push_back(new Chromosome(this->ctrl, shuffledSet, rng, false));		
		}

		if(check_interrupt()) {
			throw InterruptException();
		}
	}
	
	for(i = this->ctrl.numGenerations; i > 0; --i) {
		/*
		 *
		 * Transform the fitness map of the current generation to start at 0
		 * and copy old generation to new generation
		 */
		sumFitness = 0.0;
		IF_DEBUG(Rcout << "Fitness map: ")
		
		for(j = 0; j < this->ctrl.populationSize; ++j) {
			*(this->currentGeneration[j]) = *(newGeneration[j]);
			sumFitness += (this->currentGeneration[j]->getFitness() - minFitness);
			this->currentGenFitnessMap[j] = sumFitness;
			
			IF_DEBUG(Rcout << sumFitness << " | ")
		}
		IF_DEBUG(Rcout << std::endl)

		minFitness = 0.0;
		
		Rcpp::Rcout << "Unique chromosomes: " << this->countUniques() << std::endl;
		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}
		
		child1It = newGeneration.begin();
		child2It = newGeneration.rbegin();
		
		while(child1It != child2It.base() && child1It != (child2It + 1).base()) {
			tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(rng(0.0, sumFitness));
			tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(rng(0.0, sumFitness));

			tmpChromosome1->mateWith(*tmpChromosome2, rng, *(*child1It), *(*child2It));
			
			minParentFitness = ((tmpChromosome1->getFitness() > tmpChromosome2->getFitness()) ? tmpChromosome1->getFitness() : tmpChromosome2->getFitness());

			/*
			 * If both children have no variables, mate again
			 */
			while((*child1It)->getVariableCount() == 0 && (*child2It)->getVariableCount() == 0) {
				tmpChromosome1->mateWith(*tmpChromosome2, rng, *(*child1It), *(*child2It));
			}

			if((*child1It)->getVariableCount() == 0) {
				delete *child1It;
				*child1It = new Chromosome(**child2It);
			} else if((*child2It)->getVariableCount() == 0) {
				delete *child2It;
				*child2It = new Chromosome(**child1It);
			}

			this->evaluator.evaluate(**child1It);
			this->evaluator.evaluate(**child2It);
			// Make sure the first child is "better" than the second child
			if((*child1It)->getFitness() < (*child2It)->getFitness()) {
				std::swap(*child1It, *child2It);
			}
			
			IF_DEBUG(
				Rcout << "Mating chromosomes " << std::endl << *tmpChromosome1 << " and" << std::endl << *tmpChromosome2 << std::endl
					 << "with minimal fitness " << minParentFitness << std::endl << "First two proposals have fitness " << (*child1It)->getFitness() << " / " << (*child2It)->getFitness() << std::endl;
			)
			
			// At least the first child should be better than the worse parent
			matingTries = 0;
			while(((*child1It)->getFitness() < minParentFitness) && (++matingTries < this->ctrl.maxMatingTries)) {
				tmpChromosome1->mateWith(*tmpChromosome2, rng, *proposalChild1, *proposalChild2);
				
				/*
				 * After mating a chromosome may have no variables at all, so we need to check if the variable count is
				 * greater than 0, otherwise the evaluation step would fail
				 */
				if(proposalChild1->getVariableCount() > 0) {
					if(this->evaluator.evaluate(*proposalChild1) > (*child2It)->getFitness()) { // better as 2nd child
						if(proposalChild1->getFitness() > (*child1It)->getFitness()) { // even better as 1st child
							std::swap(*child1It, *child2It);
							delete *child1It;
							*child1It = new Chromosome(*proposalChild1);
						} else {
							delete *child2It;
							*child2It = new Chromosome(*proposalChild1);
						}
					}
				}
				
				// Check 2nd new child
				if(proposalChild2->getVariableCount() > 0) {
					if(this->evaluator.evaluate(*proposalChild2) > (*child2It)->getFitness()) { // better as 2nd child
						if(proposalChild2->getFitness() > (*child1It)->getFitness()) { // even better as 1st child
							std::swap(*child1It, *child2It);
							delete *child1It;
							*child1It = new Chromosome(*proposalChild2);
						} else {
							delete *child2It;
							*child2It = new Chromosome(*proposalChild2);
						}
					}
				}
				
				IF_DEBUG(
					Rcout << "Proposed children have fitness: " << proposalChild1->getFitness() << " / " << proposalChild2->getFitness() << std::endl
					<< "Currently selected children have fitness: " << (*child1It)->getFitness() << " / " << (*child2It)->getFitness() << std::endl;
				)
			}
			
			child1Mutated = (*child1It)->mutate(rng);
			child2Mutated = (*child2It)->mutate(rng);
			
			/*
			 * Search for identical chromosomes in the previously generated
			 * chromsomes
			 * The manual loop is faster than find_if because both children
			 * can be checked at once
			 */
			child1Dup = false;
			child2Dup = (**child1It == **child2It);
			findIt = newGeneration.begin();
			revFindIt = newGeneration.rbegin();

			while(findIt != child1It && (child1Dup == false || child2Dup == false)) {
				if(child1Dup == false && (**findIt == **child1It)) {
					child1Dup = true;
				}
				if(child2Dup == false && (**findIt == **child2It)) {
					child2Dup = true;
				}
				++findIt;
			}
			while(revFindIt != child2It && (child1Dup == false || child2Dup == false)) {
				if(child1Dup == false && (**revFindIt == **child1It)) {
					child1Dup = true;
				}
				if(child2Dup == false && (**revFindIt == **child2It)) {
					child2Dup = true;
				}
				++revFindIt;
			}
			
			if(child1Dup == false || (++child1Tries > Population::MAX_DUPLICATE_TRIES)) {
				if(child1Mutated == true) {
					this->evaluator.evaluate(**child1It);
				}

				if((*child1It)->getFitness() < minFitness) {
					minFitness = (*child1It)->getFitness();
				}

				this->addChromosomeToElite(**child1It);
				if(this->ctrl.verbosity >= MORE_VERBOSE) {
					this->printChromosomeFitness(Rcout, **child1It);
				}

				++child1It;

				IF_DEBUG(
					if(child1Tries > MAX_DUPLICATE_TRIES) {
						Rcpp::Rcout << "Needed " << (int) child1Tries << " tries to find unique chromosome" << std::endl;
					}
				)
				child1Tries = 0;
			}
			
			if(child2Dup == false || (++child2Tries > Population::MAX_DUPLICATE_TRIES)) {
				if(child2Mutated == true) {
					this->evaluator.evaluate(**child2It);
				}

				if((*child2It)->getFitness() < minFitness) {
					minFitness = (*child2It)->getFitness();
				}

				this->addChromosomeToElite(**child2It);
				if(this->ctrl.verbosity >= MORE_VERBOSE) {
					this->printChromosomeFitness(Rcout, **child2It);
				}

				++child2It;

				IF_DEBUG(
					if(child2Tries > MAX_DUPLICATE_TRIES) {
						Rcpp::Rcout << "Needed " << (int) child2Tries << " tries to find unique chromosome" << std::endl;
					}
				)
				child2Tries = 0;
			}
			
			if(check_interrupt()) {
				throw InterruptException();
			}
		}
	}
	
	delete proposalChild1;
	delete proposalChild2;
	for(std::vector<Chromosome*>::iterator it = newGeneration.begin(); it != newGeneration.end(); ++it) {
		delete *it;
	}
}