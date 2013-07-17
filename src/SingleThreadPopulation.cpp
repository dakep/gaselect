//
//  SingleThreadPopulation.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//
#include "config.h"

#include <vector>
#include <algorithm>
#include <RcppArmadillo.h>

#include "UnifGenerator_0_1.h"
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

SingleThreadPopulation::SingleThreadPopulation(const Control &ctrl, ::Evaluator &evaluator) : Population(ctrl, evaluator) {}


void SingleThreadPopulation::run() {
	int i = 0, j = 0;
	int popSizeHalf = this->ctrl.populationSize / 2;
	VariablePositionPopulation varPosPop(this->ctrl.chromosomeSize);
	
	double sumFitness = 0.0;
	double minFitness = 0.0;
	
	//	Chromosome tmpChromosome1, tmpChromosome2;
	std::vector<double> newFitnessMap;
	std::vector<Chromosome> newGeneration;
	std::vector<double>::iterator fitnessMapIt;
	UnifGenerator_0_1 unifGen;
	
	newGeneration.reserve(this->ctrl.populationSize);
	newFitnessMap.reserve(this->ctrl.populationSize);
	
	if(this->ctrl.verbosity > OFF) {
		Rcout << "Generating initial population" << std::endl;
	}
	
	for(i = this->ctrl.populationSize; i > 0; --i) {
		Chromosome tmpChromosome1(this->ctrl, varPosPop, unifGen);
		
		this->currentGenFitnessMap.push_back(this->evaluator.evaluate(tmpChromosome1));
		if(tmpChromosome1.getFitness() < minFitness) {
			minFitness = tmpChromosome1.getFitness();
		}
		
		if(this->ctrl.verbosity >= MORE_VERBOSE) {
			this->printChromosomeFitness(Rcout, tmpChromosome1);
		}
		this->addChromosomeToElite(tmpChromosome1);
		this->currentGeneration.push_back(tmpChromosome1);
		
		if(check_interrupt()) {
			throw InterruptException();
		}
	}
		
	for(i = this->ctrl.numGenerations; i > 0; --i) {
		/*
		 *
		 * Transform the fitness map of the current generation to start at 0
		 *
		 */
		sumFitness = 0.0;
		std::vector<double>::iterator fitnessMapIt;
		IF_DEBUG(Rcout << "Fitness map: ")
		
		for(fitnessMapIt = this->currentGenFitnessMap.begin(); fitnessMapIt != this->currentGenFitnessMap.end(); ++fitnessMapIt) {
			sumFitness += (*fitnessMapIt) - minFitness;
			(*fitnessMapIt) = sumFitness;
			
			IF_DEBUG(Rcout << sumFitness << " | ")
		}
		IF_DEBUG(Rcout << std::endl)

		minFitness = 0.0;
		
		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}
		
		// Copulate and mutate
		for(j = 0; j < popSizeHalf; ++j) {
			Chromosome tmpChromosome1 = this->drawChromosomeFromCurrentGeneration(unifGen() * sumFitness);
			Chromosome tmpChromosome2 = this->drawChromosomeFromCurrentGeneration(unifGen() * sumFitness);
			std::vector<Chromosome> children = tmpChromosome1.mateWith(tmpChromosome2, unifGen);
			uint8_t matingTries = 0;
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

			this->evaluator.evaluate(children[0]);
			this->evaluator.evaluate(children[1]);
			// Make sure the first child is "better" than the second child
			if(children[0].getFitness() < children[1].getFitness()) {
				std::swap(children[1], children[0]);
			}
			
			
			IF_DEBUG(
				Rcout << "Mating chromosomes " << std::endl << tmpChromosome1 << " and" << std::endl << tmpChromosome2 << std::endl
					 << "with minimal fitness " << minParentFitness << std::endl << "First two proposals have fitness " << children[0].getFitness() << " / " << children[1].getFitness() << std::endl;
			)
			
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
			
			newFitnessMap.push_back(children[0].getFitness());
			newFitnessMap.push_back(children[1].getFitness());
			
			if(children[0].getFitness() < minFitness) {
				minFitness = children[0].getFitness();
			}
			
			if(children[1].getFitness() < minFitness) {
				minFitness = children[1].getFitness();
			}
			
			if(this->ctrl.verbosity >= MORE_VERBOSE) {
				this->printChromosomeFitness(Rcout, children[0]);
				this->printChromosomeFitness(Rcout, children[1]);
			}
			
			this->addChromosomeToElite(children[0]);
			this->addChromosomeToElite(children[1]);
			
			newGeneration.push_back(children[0]);
			newGeneration.push_back(children[1]);
			
			if(check_interrupt()) {
				throw InterruptException();
			}
		}
		
		// Housekeeping
		// first delete old chromsomes
		//		this->cleanCurrentGeneration();
		this->currentGeneration = newGeneration;
		newGeneration.clear();
		
		this->currentGenFitnessMap = newFitnessMap;
		newFitnessMap.clear();
	}
}