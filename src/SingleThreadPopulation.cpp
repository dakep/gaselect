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
#define IF_DEBUG(expr) if(this->ctrl.verbosity == DEBUG_GA || this->ctrl.verbosity == DEBUG_ALL) { expr; }
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
	uint8_t matingTries = 0;
	ShuffledSet shuffledSet(this->ctrl.chromosomeSize);
	RNG rng(this->seed);
	
	double sumFitness = 0.0;
	double minFitness = 0.0;
	double minParentFitness = 0.0;
	
	ChromosomeVec newGeneration;
	
	Chromosome* tmpChromosome1;
	Chromosome* tmpChromosome2;
	ChromosomeVecIter child1It;
	ChromosomeVecRevIter child2It;
	Chromosome* proposalChild1 = new Chromosome(this->ctrl, shuffledSet, rng, false);
	Chromosome* proposalChild2 = new Chromosome(this->ctrl, shuffledSet, rng, false);
	uint8_t child1Tries = 0;
	uint8_t child2Tries = 0;
	bool child1Mutated, child2Mutated;
	std::pair<bool, bool> duplicated;
	
	newGeneration.reserve(this->ctrl.populationSize);
	
	if(this->ctrl.verbosity > OFF) {
		Rcout << "Generating initial population" << std::endl;
	}
	
	while(newGeneration.size() < this->ctrl.populationSize) {
		tmpChromosome1 = new Chromosome(this->ctrl, shuffledSet, rng);
		
		/* Check if chromosome is already in the initial population */
		if(std::find_if(newGeneration.begin(), newGeneration.end(), CompChromsomePtr(tmpChromosome1)) == newGeneration.end()) {
			this->evaluator.evaluate(*tmpChromosome1);

			if(tmpChromosome1->getFitness() < minFitness) {
				minFitness = tmpChromosome1->getFitness();
			}
			
			this->addChromosomeToElite(*tmpChromosome1);
			
			newGeneration.push_back(tmpChromosome1);
		} else {
			delete tmpChromosome1;
		}

		if(check_interrupt()) {
			throw InterruptException();
		}
	}

	this->initCurrentGeneration(shuffledSet, rng);

	/*
	 * Transform the fitness map of the current generation to start at 0
	 * and copy old generation to new generation
	 */
	sumFitness = this->updateCurrentGeneration(newGeneration, minFitness);

	if(this->ctrl.verbosity >= VERBOSE && this->ctrl.verbosity != DEBUG_EVAL) {
		this->printCurrentGeneration();
	}
	
	for(i = this->ctrl.numGenerations; i > 0; --i) {
		minFitness = 0.0;

		IF_DEBUG(
			Rcpp::Rcout << "Unique chromosomes: " << this->countUniques() << std::endl;
		)

		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}
		
		child1It = newGeneration.begin();
		child2It = newGeneration.rbegin();
		
		while(child1It != child2It.base() && child1It != (child2It + 1).base()) {
			tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(rng(0.0, sumFitness));
			do {
				tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(rng(0.0, sumFitness));
			} while(tmpChromosome1 == tmpChromosome2);

			tmpChromosome1->mateWith(*tmpChromosome2, rng, *(*child1It), *(*child2It));
			
			minParentFitness = ((tmpChromosome1->getFitness() > tmpChromosome2->getFitness()) ? tmpChromosome1->getFitness() : tmpChromosome2->getFitness());

			/*
			 * If any of the two children has no variables, mate again
			 */
			while((*child1It)->getVariableCount() == 0 || (*child2It)->getVariableCount() == 0) {
				tmpChromosome1->mateWith(*tmpChromosome2, rng, *(*child1It), *(*child2It));
			}

			this->evaluator.evaluate(**child1It);
			this->evaluator.evaluate(**child2It);

			// Make sure the first child is "better" than the second child
			if((*child1It)->getFitness() < (*child2It)->getFitness()) {
				std::swap(*child1It, *child2It);
			}

			matingTries = 0;
			while(((*child1It)->getFitness() < minParentFitness) && (matingTries++ < this->ctrl.maxMatingTries)) {
				tmpChromosome1->mateWith(*tmpChromosome2, rng, *proposalChild1, *proposalChild2);
				
				if((*child1It)->getFitness() < minParentFitness && matingTries++ < this->ctrl.maxMatingTries) {
					/*
					 * After mating a chromosome may have no variables at all, so we need to check if the variable count is
					 * greater than 0, otherwise the evaluation step would fail
					 * The better child is worse than the worst parent and we have some tries left
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
			}
			
			if((*child1It)->getFitness() < (minParentFitness - this->ctrl.badSolutionThreshold * fabs(minParentFitness))) {
				/*
				 * The fitness of the better child is more than x% less than the worst parent's
				 * fitness so cancel mating of the two parents and choose two new parents
				 */
				continue;
			}

			child1Mutated = (*child1It)->mutate(rng);
			child2Mutated = (*child2It)->mutate(rng);
			
			/*
			 * Search for identical chromosomes in the previously generated
			 * chromsomes
			 * The manual loop is faster than find_if because both children
			 * can be checked at once
			 */
			duplicated = Population::checkDuplicated(newGeneration.begin(), newGeneration.rbegin(), child1It, child2It);
			
			if(duplicated.first == false || (++child1Tries > this->ctrl.maxDuplicateEliminationTries)) {
				if(child1Mutated == true) {
					this->evaluator.evaluate(**child1It);
				}

				if((*child1It)->getFitness() < minFitness) {
					minFitness = (*child1It)->getFitness();
				}

				this->addChromosomeToElite(**child1It);
				if(this->ctrl.verbosity >= VERBOSE) {
					this->printChromosomeFitness(Rcout, **child1It);
				}

				++child1It;

				IF_DEBUG(
					if(child1Tries > 0) {
						Rcpp::Rcout << "Needed " << (int) child1Tries << " tries to find unique chromosome" << std::endl;
					}
				)
				child1Tries = 0;
			}
			
			if(duplicated.second == false || (++child2Tries > this->ctrl.maxDuplicateEliminationTries)) {
				if(child2Mutated == true) {
					this->evaluator.evaluate(**child2It);
				}

				if((*child2It)->getFitness() < minFitness) {
					minFitness = (*child2It)->getFitness();
				}

				this->addChromosomeToElite(**child2It);
				if(this->ctrl.verbosity >= VERBOSE) {
					this->printChromosomeFitness(Rcout, **child2It);
				}

				++child2It;

				IF_DEBUG(
					if(child2Tries > 0) {
						Rcpp::Rcout << "Needed " << (int) child2Tries << " tries to find unique chromosome" << std::endl;
					}
				)
				child2Tries = 0;
			}
			
			if(check_interrupt()) {
				throw InterruptException();
			}
		}

		/*
		 * Transform the fitness map of the current generation to start at 0
		 * and copy old generation to new generation
		 */		
		sumFitness = this->updateCurrentGeneration(newGeneration, minFitness, false);

		if(this->ctrl.verbosity >= VERBOSE && this->ctrl.verbosity != DEBUG_EVAL) {
			this->printCurrentGeneration();
		}

	}
	delete proposalChild1;
	delete proposalChild2;
	for(ChromosomeVecIter it = newGeneration.begin(); it != newGeneration.end(); ++it) {
		delete *it;
	}
}