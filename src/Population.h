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

#include "RNG.h"
#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->ctrl.verbosity >= DEBUG_VERBOSE) { expr; }
#else
#define IF_DEBUG(expr)
#endif

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

	Population(const Control &ctrl, ::Evaluator &evaluator, const std::vector<uint32_t> &seed) : ctrl(ctrl), evaluator(evaluator), seed(seed) {
		this->currentGeneration.reserve(this->ctrl.populationSize);
		this->currentGenFitnessMap.reserve(this->ctrl.populationSize);
		this->minEliteFitness = 0.0;
		this->cutoffPoint = this->ctrl.populationSize * this->ctrl.cutoffQuantile / 100;
	}

	virtual ~Population() {
		for(std::vector<Chromosome*>::iterator it = this->currentGeneration.begin(); it != this->currentGeneration.end(); ++it) {
			delete *it;
		}
	}

	virtual void run() = 0;

	std::multiset<Chromosome, Population::ChromosomeComparator> getResult() const {
		std::multiset<Chromosome, Population::ChromosomeComparator> result;
		
		for(std::vector<Chromosome*>::const_iterator it = this->currentGeneration.begin(); it != this->currentGeneration.end(); ++it) {
			result.insert(**it);
		}
		
		if(this->ctrl.elitism > 0 && this->elite.size() > 0) {
			result.insert(this->elite.begin(), this->elite.end());
		}
		
		return result;
	}
private:
	uint16_t cutoffPoint;
protected:
	const Control& ctrl;
	::Evaluator& evaluator;
	const std::vector<uint32_t> &seed;

	std::multiset<Chromosome, Population::ChromosomeComparator> elite;
	std::vector<Chromosome*> currentGeneration;
	std::vector<double> currentGenFitnessMap;
	double minEliteFitness;
	
	/**
	 * Pick a chromosome from the current generation at random
	 * where the probability to pick a chromosome is taken from
	 * the currentGenFitnessMap
	 */
	Chromosome* drawChromosomeFromCurrentGeneration(double rand) {
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
	
	static bool comp(Chromosome* c1, Chromosome* c2) {
		return ((*c1) == (*c2));
	}
	
	uint16_t countUniques() {
		std::vector<Chromosome*> gen = this->currentGeneration;
		std::vector<Chromosome*>::iterator end = std::unique(gen.begin(), gen.end(), Population::comp);
		return std::distance(gen.begin(), end);
	};
	
	double getQuantileFitness() {
		std::vector<double> copyFitnessMap = this->currentGenFitnessMap;
		std::vector<double>::iterator cutoffIt = copyFitnessMap.begin() + this->cutoffPoint;
		std::nth_element(copyFitnessMap.begin(), cutoffIt, copyFitnessMap.end());
		return (*cutoffIt);
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
			
			IF_DEBUG(
				Rcpp::Rcout << "Adding chromosome to elite. New minimum fitness for elite is " << this->minEliteFitness << std::endl;
			)
		}
	};
};

typedef std::multiset<Chromosome, Population::ChromosomeComparator> SortedChromosomes;
#endif
