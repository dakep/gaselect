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

Population::Population(const Control &ctrl, const ::Evaluator &evaluator) : ctrl(ctrl), evaluator(&evaluator) {
	// Initialize a vector of doubles that is used to
	// map a 0-1 uniform random variable to the appropriate
	// chromosome.
	this->fitnessMap.reserve(this->ctrl.populationSize);

	// initialize original population (generation 0) totally randomly
	this->currentGeneration.reserve(this->ctrl.populationSize);

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

Chromosome& Population::getChromosomeFromFitnessMap(double rand) {
	int imin = 0, imax = this->ctrl.populationSize;
	int imid = 0;
	bool found = false;

	while(imax > imin && !found) {
		imid = imin + ((imax - imin) / 2);

		if(this->fitnessMap[imid] < rand) {
			imin = imid + 1;
		} else if (this->fitnessMap[imid] >= rand) {
			if(imid > 0 && this->fitnessMap[imid - 1] >= rand) {
				imax = imid - 1;
			} else {
				imax = imin; // break loop
			}
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

void Population::run() {
	int i = 0, j = 0;
	int popSizeHalf = this->ctrl.populationSize / 2;
	VariablePositionPopulation varPosPop(this->ctrl.chromosomeSize);

	double sumFitness = 0.0;
	double minFitness = 0.0;

//	Chromosome tmpChromosome1, tmpChromosome2;
	std::vector<double> newFitnessMap;
	std::vector<Chromosome> newGeneration;
	std::vector<double>::iterator fitnessMapIt;
	Rcpp::stats::UnifGenerator__0__1 unifGen;

	newGeneration.reserve(this->ctrl.populationSize);
	newFitnessMap.reserve(this->ctrl.populationSize);

	if(this->ctrl.verbosity > OFF) {
		Rcout << "Generating initial population:" << std::endl;
	}

	for(i = this->ctrl.populationSize; i > 0; --i) {
		Chromosome tmpChromosome1(this->ctrl, varPosPop);

		this->fitnessMap.push_back(this->evaluator->evaluate(tmpChromosome1));
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

#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << "Fitness map: ";
	}
#endif

	for(fitnessMapIt = this->fitnessMap.begin(); fitnessMapIt != this->fitnessMap.end(); ++fitnessMapIt) {
		sumFitness += (*fitnessMapIt) - minFitness;
		(*fitnessMapIt) = sumFitness;

#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.verbosity == DEBUG_VERBOSE) {
			Rcout << sumFitness << " | ";
		}
#endif
	}
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << std::endl;
	}
#endif

	for(i = this->ctrl.numGenerations; i > 0; --i) {
		if(this->ctrl.verbosity > OFF) {
			Rcout << "Generating generation " << (this->ctrl.numGenerations - i + 1) << std::endl;
		}

		minFitness = 0.0;

		// Copulate and mutate
		for(j = 0; j < popSizeHalf; ++j) {
			Chromosome tmpChromosome1 = this->getChromosomeFromFitnessMap(unifGen() * sumFitness);
			Chromosome tmpChromosome2 = this->getChromosomeFromFitnessMap(unifGen() * sumFitness);
			std::vector<Chromosome> children = tmpChromosome1.mateWith(tmpChromosome2);
			uint8_t cnt = 0;
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
				<< "First two proposals have fitness " << children[0].getFitness() << " / " << children[2].getFitness() << std::endl;
			}
#endif
			
			// At least the first child should be better than the worse parent
			while((children[0].getFitness() < minParentFitness) && (++cnt < this->ctrl.maxMatingTries)) {
				std::vector<Chromosome> proposalChildren = tmpChromosome1.mateWith(tmpChromosome2);
				if(this->evaluator->evaluate(proposalChildren[0]) > children[1].getFitness()) { // better as 2nd child
					if(proposalChildren[0].getFitness() > children[0].getFitness()) { // even better as 1st child
						children[1] = children[0];
						children[0] = proposalChildren[0];
					} else {
						children[1] = proposalChildren[0];
					}
				}

				// Check 2nd new child
				if(this->evaluator->evaluate(proposalChildren[1]) > children[1].getFitness()) { // better as 2nd child
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
			
			children[0].mutate();
			children[1].mutate();

			newFitnessMap.push_back(this->evaluator->evaluate(children[0]));
			newFitnessMap.push_back(this->evaluator->evaluate(children[1]));

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

#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.verbosity == DEBUG_VERBOSE) {
			Rcout << "Fitness map: ";
		}
#endif

		sumFitness = 0.0;

		for(fitnessMapIt = newFitnessMap.begin(); fitnessMapIt != newFitnessMap.end(); ++fitnessMapIt) {
			sumFitness += (*fitnessMapIt) - minFitness;
			(*fitnessMapIt) = sumFitness;

#ifdef ENABLE_DEBUG_VERBOSITY
			if(this->ctrl.verbosity == DEBUG_VERBOSE) {
				Rcout << sumFitness << " | ";
			}
#endif
		}
#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.verbosity == DEBUG_VERBOSE) {
			Rcout << std::endl;
		}
#endif

		// Housekeeping
		// first delete old chromsomes
//		this->cleanCurrentGeneration();
		this->currentGeneration = newGeneration;
		newGeneration.clear();

		this->fitnessMap = newFitnessMap;
		newFitnessMap.clear();
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


//inline void Population::cleanCurrentGeneration() {
//	std::vector<Chromosome>::iterator rmGenerationIter = this->currentGeneration.begin();
//
//	// Delete all chromosomes from the last generation from stack
//	for(; rmGenerationIter != this->currentGeneration.end(); ++rmGenerationIter) {
//		delete (*rmGenerationIter);
//	}
//}
