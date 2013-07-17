//
//  Population.h
//
//

#ifndef GenAlgPLS_Population_h
#define GenAlgPLS_Population_h

#include "config.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <set>

#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"

class Population {
public:
	struct ChromosomeComparator {
		bool operator() (const Chromosome &lhs, const Chromosome &rhs) const {
			return rhs.isFitterThan(lhs);
		}
	};
	
	class InterruptException : public std::exception {
	public:
		InterruptException() {};
		virtual ~InterruptException() throw() {};
		virtual const char* what() const throw() {
			return "Interrupted";
		}
	};

	Population(const Control &ctrl, ::Evaluator &evaluator) : ctrl(ctrl), evaluator(evaluator) {
		this->currentGeneration.reserve(this->ctrl.populationSize);
		this->currentGenFitnessMap.reserve(this->ctrl.populationSize);
		this->minEliteFitness = 0.0;
	}
	virtual ~Population() {};

	virtual void run() = 0;

	std::multiset<Chromosome, Population::ChromosomeComparator> getResult() const {
		std::multiset<Chromosome, Population::ChromosomeComparator> result(this->currentGeneration.begin(), this->currentGeneration.end());
		
		if(this->ctrl.elitism > 0 && this->elite.size() > 0) {
			result.insert(this->elite.begin(), this->elite.end());
		}
		
		return result;
	}
protected:
	const Control& ctrl;
	::Evaluator& evaluator;

	std::multiset<Chromosome, Population::ChromosomeComparator> elite;
	std::vector<Chromosome> currentGeneration;
	std::vector<double> currentGenFitnessMap;
	double minEliteFitness;
	
	/**
	 * Pick a chromosome from the current generation at random
	 * where the probability to pick a chromosome is taken from
	 * the currentGenFitnessMap
	 */
	Chromosome &drawChromosomeFromCurrentGeneration(double rand) {
		int imin = 0, imax = this->ctrl.populationSize - 1;
		int imid = 0;
		
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
	};

	std::ostream& printChromosomeFitness(std::ostream &os, Chromosome &ch) {
		os << ch << TAB_DELIMITER << std::fixed
		<< std::setw(WIDTH) << std::setprecision(PRECISION) << ch.getFitness() << std::endl;
		
		return os;
	};
	
	void addChromosomeToElite(Chromosome &ch) {
		// Add chromosome to elite if better than the worst elite-chromosome
		if(this->ctrl.elitism > 0 && (ch.getFitness() > this->minEliteFitness || this->elite.size() < this->ctrl.elitism)) {
			if(this->elite.size() >= this->ctrl.elitism) {
				this->elite.erase(this->elite.begin());
			}

			this->elite.insert(ch); // The insert copies the chromosome!
			
			this->minEliteFitness = this->elite.begin()->getFitness();
			
			if(this->ctrl.verbosity >= MORE_VERBOSE) {
				Rcpp::Rcout << "Adding chromosome to elite. New minimum fitness for elite is " << this->minEliteFitness << std::endl;
			}
		}
	};
};

typedef std::multiset<Chromosome, Population::ChromosomeComparator> SortedChromosomes;
#endif
