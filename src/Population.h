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
#include <utility>

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
protected:
	typedef std::vector<Chromosome*> ChromosomeVec;
	typedef std::vector<Chromosome*>::iterator ChromosomeVecIter;
	typedef std::vector<Chromosome*>::reverse_iterator ChromosomeVecRevIter;

public:
	class ChromosomeComparator : public std::binary_function<Chromosome, Chromosome, bool> {
	public:
		bool operator() (const Chromosome &lhs, const Chromosome &rhs) const {
			return rhs.isFitterThan(lhs);
		}
	};
	
	typedef std::set<Chromosome, Population::ChromosomeComparator> SortedChromosomes;
	
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
	}
	
	virtual ~Population() {
		for(ChromosomeVecIter it = this->currentGeneration.begin(); it != this->currentGeneration.end(); ++it) {
			delete *it;
		}
	}
	
	virtual void run() = 0;
	
	inline SortedChromosomes getResult() const {
		SortedChromosomes result(this->elite);
		
		/*
		 * result.insert(begin(), end()) can not be used because we need to
		 * insert the VALUE and not the pointer to a chromosome
		 */
		for(ChromosomeVec::const_iterator it = this->currentGeneration.begin(); it != this->currentGeneration.end(); ++it) {
			result.insert(**it);
		}
		
		return result;
	}
protected:
	const Control& ctrl;
	::Evaluator& evaluator;
	const std::vector<uint32_t> &seed;
	
	SortedChromosomes elite;
	ChromosomeVec currentGeneration;
	std::vector<double> currentGenFitnessMap;
	double minEliteFitness;
	
	/**
	 * Pick a chromosome from the current generation at random
	 * where the probability to pick a chromosome is taken from
	 * the currentGenFitnessMap
	 */
	inline Chromosome* drawChromosomeFromCurrentGeneration(double rand) {
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

#ifdef ENABLE_DEBUG_VERBOSITY
	static bool compEqual(Chromosome* c1, Chromosome* c2) {
		return ((*c1) == (*c2));
	}
	
	static bool compLT(const Chromosome* const c1, const Chromosome* const c2) {
		return c1->isFitterThan(*c2);
	}
	
	inline uint16_t countUniques() {
		std::vector<Chromosome*> gen = this->currentGeneration;
		std::sort(gen.begin(), gen.end(), Population::compLT);
		return std::distance(gen.begin(), std::unique(gen.begin(), gen.end(), Population::compEqual));
	};
#endif
	
	inline std::ostream& printChromosomeFitness(std::ostream &os, Chromosome &ch) {
		os << ch << TAB_DELIMITER << std::fixed
		<< std::setw(WIDTH) << std::setprecision(PRECISION) << ch.getFitness() << std::endl;
		
		return os;
	};
	
	inline void addChromosomeToElite(Chromosome &ch) {
		/* Add chromosome to elite if better than the worst elite-chromosome */
		if(this->ctrl.elitism > 0 && (ch.getFitness() > this->minEliteFitness || this->elite.size() < this->ctrl.elitism)) {
			/* Insert a copy of the chromosome. If the chromosome is a duplicate, it is not inserted */
			this->elite.insert(ch);
			
			/*
			 * If the chromosome was inserted and the elite-set was already `full`
			 * the last worst chromosome is removed
			 */
			if(this->elite.size() > this->ctrl.elitism) {
				this->elite.erase(this->elite.begin());
			}
			
			this->minEliteFitness = this->elite.begin()->getFitness();
			
			IF_DEBUG(
				Rcpp::Rcout << "Adding chromosome to elite. New minimum fitness for elite is " << this->minEliteFitness << std::endl;
			)
		}
	};
	
	inline static std::pair<bool, bool> checkDuplicated(ChromosomeVecIter begin, ChromosomeVecRevIter rbegin, const ChromosomeVecIter &child1It, const ChromosomeVecRevIter &child2It) {
		std::pair<bool, bool> duplicated(false, **child1It == **child2It);
		
		while(begin != child1It && (duplicated.first == false || duplicated.second == false)) {
			if(duplicated.first == false && (**begin == **child1It)) {
				duplicated.first = true;
			}
			if(duplicated.second == false && (**begin == **child2It)) {
				duplicated.second = true;
			}
			++begin;
		}
		while(rbegin != child2It && (duplicated.first == false || duplicated.second == false)) {
			if(duplicated.first == false && (**rbegin == **child1It)) {
				duplicated.first = true;
			}
			if(duplicated.second == false && (**rbegin == **child2It)) {
				duplicated.second = true;
			}
			++rbegin;
		}
		
		return duplicated;
	};

	class CompChromsomePtr : public std::unary_function<Chromosome*, bool> {
	public:
		CompChromsomePtr(const Chromosome* const ch) : ch(ch) {}
		bool operator()(const Chromosome* const ch) const {
			return ((*this->ch) == (*ch));
		}
		
	private:
		const Chromosome* const ch;
	};
};
#endif
