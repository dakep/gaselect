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

SingleThreadPopulation::SingleThreadPopulation(const Control &ctrl, ::Evaluator &evaluator, RNG &rng) : Population(ctrl, evaluator, rng) {}


void SingleThreadPopulation::run() {
	int i = 0, j = 0;
	uint8_t matingTries = 0;
	int popSizeHalf = this->ctrl.populationSize / 2;
	VariablePositionPopulation varPosPop(this->ctrl.chromosomeSize);
	
	double sumFitness = 0.0;
	double minFitness = 0.0;
	double minParentFitness = 0.0;
	
	//	Chromosome tmpChromosome1, tmpChromosome2;
	std::vector<double> newFitnessMap;
	std::vector<Chromosome*> newGeneration;
	std::vector<double>::iterator fitnessMapIt;
	
	Chromosome* tmpChromosome1;
	Chromosome* tmpChromosome2;
	Chromosome** child1;
	Chromosome** child2;
	Chromosome* proposalChild1 = new Chromosome(this->ctrl, varPosPop, this->rng, false);
	Chromosome* proposalChild2 = new Chromosome(this->ctrl, varPosPop, this->rng, false);
	
	newGeneration.reserve(this->ctrl.populationSize);
	newFitnessMap.reserve(this->ctrl.populationSize);
	
	if(this->ctrl.verbosity > OFF) {
		Rcout << "Generating initial population" << std::endl;
	}
	
	for(i = this->ctrl.populationSize; i > 0; --i) {
		tmpChromosome1 = new Chromosome(this->ctrl, varPosPop, this->rng);
		
		this->currentGenFitnessMap.push_back(this->evaluator.evaluate(*tmpChromosome1));
		if(tmpChromosome1->getFitness() < minFitness) {
			minFitness = tmpChromosome1->getFitness();
		}
		
		if(this->ctrl.verbosity >= MORE_VERBOSE) {
			this->printChromosomeFitness(Rcout, *tmpChromosome1);
		}
		this->addChromosomeToElite(*tmpChromosome1);
		this->currentGeneration.push_back(tmpChromosome1);
		newGeneration.push_back(new Chromosome(this->ctrl, varPosPop, this->rng, false));
		
		if(check_interrupt()) {
			throw InterruptException();
		}
	}
	
	double scale = 0.0;
		
	for(i = this->ctrl.numGenerations; i > 0; --i) {
		/*
		 *
		 * Transform the fitness map of the current generation to start at 0
		 *
		 */
		sumFitness = 0.0;
		IF_DEBUG(Rcout << "Fitness map: ")
		
		scale = (this->ctrl.cutoffQuantile > 0 ? this->getQuantileFitness() : minFitness);

		for(fitnessMapIt = this->currentGenFitnessMap.begin(); fitnessMapIt != this->currentGenFitnessMap.end(); ++fitnessMapIt) {
			if((*fitnessMapIt) < scale) {
				(*fitnessMapIt) = 0.0;
			} else {
				sumFitness += (*fitnessMapIt) - scale;
				(*fitnessMapIt) = sumFitness;
			}
//			sumFitness += (*fitnessMapIt) - minFitness;
//			(*fitnessMapIt) = sumFitness;
			
			IF_DEBUG(Rcout << sumFitness << " | ")
		}
		IF_DEBUG(Rcout << std::endl)

		minFitness = 0.0;
		
		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}
		
		for(j = 0; j < popSizeHalf; ++j) {
			tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(this->rng(0.0, sumFitness));
			tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(this->rng(0.0, sumFitness));
			
			child1 = &(newGeneration[2 * j]);
			child2 = &(newGeneration[(2 * j) + 1]);

			tmpChromosome1->mateWith(*tmpChromosome2, this->rng, **child1, **child2);
			
			minParentFitness = ((tmpChromosome1->getFitness() > tmpChromosome2->getFitness()) ? tmpChromosome1->getFitness() : tmpChromosome2->getFitness());

			/*
			 * If both children have no variables, mate again
			 */
			while((*child1)->getVariableCount() == 0 && (*child2)->getVariableCount() == 0) {
				tmpChromosome1->mateWith(*tmpChromosome2, this->rng, **child1, **child2);
			}

			if((*child1)->getVariableCount() == 0) {
				delete *child1;
				*child1 = new Chromosome(**child2);
			} else if((*child2)->getVariableCount() == 0) {
				delete *child2;
				*child2 = new Chromosome(**child1);
			}

			this->evaluator.evaluate(**child1);
			this->evaluator.evaluate(**child2);
			// Make sure the first child is "better" than the second child
			if((*child1)->getFitness() < (*child2)->getFitness()) {
				std::swap(child1, child2);
			}
			
			IF_DEBUG(
				Rcout << "Mating chromosomes " << std::endl << *tmpChromosome1 << " and" << std::endl << *tmpChromosome2 << std::endl
					 << "with minimal fitness " << minParentFitness << std::endl << "First two proposals have fitness " << (*child1)->getFitness() << " / " << (*child2)->getFitness() << std::endl;
			)
			
			// At least the first child should be better than the worse parent
			matingTries = 0;
			while(((*child1)->getFitness() < minParentFitness) && (++matingTries < this->ctrl.maxMatingTries)) {
				tmpChromosome1->mateWith(*tmpChromosome2, this->rng, *proposalChild1, *proposalChild2);
				
				/*
				 * After mating a chromosome may have no variables at all, so we need to check if the variable count is
				 * greater than 0, otherwise the evaluation step would fail
				 */
				if(proposalChild1->getVariableCount() > 0) {
					if(this->evaluator.evaluate(*proposalChild1) > (*child2)->getFitness()) { // better as 2nd child
						if(proposalChild1->getFitness() > (*child1)->getFitness()) { // even better as 1st child
							std::swap(child1, child2);
							delete *child1;
							*child1 = new Chromosome(*proposalChild1);
						} else {
							delete *child2;
							*child2 = new Chromosome(*proposalChild1);
						}
					}
				}
				
				// Check 2nd new child
				if(proposalChild2->getVariableCount() > 0) {
					if(this->evaluator.evaluate(*proposalChild2) > (*child2)->getFitness()) { // better as 2nd child
						if(proposalChild2->getFitness() > (*child1)->getFitness()) { // even better as 1st child
							std::swap(child1, child2);
							delete *child1;
							*child1 = new Chromosome(*proposalChild2);
						} else {
							delete *child2;
							*child2 = new Chromosome(*proposalChild2);
						}
					}
				}
				
				IF_DEBUG(
					Rcout << "Proposed children have fitness: " << proposalChild1->getFitness() << " / " << proposalChild2->getFitness() << std::endl
					<< "Currently selected children have fitness: " << (*child1)->getFitness() << " / " << (*child2)->getFitness() << std::endl;
				)
			}
			
			if((*child1)->mutate(this->rng) == true) {
				this->evaluator.evaluate(**child1);
			}
			if((*child2)->mutate(this->rng) == true) {
				this->evaluator.evaluate(**child2);
			}
			
			newFitnessMap.push_back((*child1)->getFitness());
			newFitnessMap.push_back((*child2)->getFitness());
			
			if((*child1)->getFitness() < minFitness) {
				minFitness = (*child1)->getFitness();
			}
			
			if((*child2)->getFitness() < minFitness) {
				minFitness = (*child2)->getFitness();
			}
			
			if(this->ctrl.verbosity >= MORE_VERBOSE) {
				this->printChromosomeFitness(Rcout, **child1);
				this->printChromosomeFitness(Rcout, **child2);
			}
			
			this->addChromosomeToElite(**child1);
			this->addChromosomeToElite(**child2);
			
			if(check_interrupt()) {
				throw InterruptException();
			}
		}
		
		/*
		 * Copy the objects from the new generation to the old generation
		 */
		for(j = this->ctrl.populationSize - 1; j >= 0; --j) {
			*(this->currentGeneration[j]) = *(newGeneration[j]);
		}
		
		this->currentGenFitnessMap = newFitnessMap;
		newFitnessMap.clear();
	}
	
	delete proposalChild1;
	delete proposalChild2;
	for(std::vector<Chromosome*>::iterator it = newGeneration.begin(); it != newGeneration.end(); ++it) {
		delete *it;
	}
}