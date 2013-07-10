//
//  Population.h
//
//

#ifndef GenAlgPLS_Population_h
#define GenAlgPLS_Population_h

#include "config.h"

#include <iostream>
#include <vector>
#include <set>

#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"

class Population {

public:
	Population(const Control &ctrl, const ::Evaluator &evaluator);
	~Population();

	void run();

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

	/**
	 * Returns the last generation and the elite (if any)
	 * Invalid once the Population object is destroyed!
	 */
	std::multiset<Chromosome, Population::ChromosomeComparator> getResult() const;

private:
	std::multiset<Chromosome, Population::ChromosomeComparator> elite;
	std::vector<Chromosome> currentGeneration;
	std::vector<Chromosome> nextGeneration;
	std::vector<double> currentGenFitnessMap;
	std::vector<double> nextGenFitnessMap;
	double sumCurrentGenFitness;
	double minCurrentGenFitness; // Minimum fitness value in the current generation
	double minEliteFitness;
	Rcpp::stats::UnifGenerator__0__1 unifGen;

	/**
	 * Pick a chromosome from the current generation at random
	 * where the probability to pick a chromosome is taken from
	 * the currentGenFitnessMap
	 */
	Chromosome &drawChromosomeFromCurrentGeneration();
	void addChromosomeToElite(Chromosome &ch);

	std::ostream& printChromosomeFitness(std::ostream &os, Chromosome &ch);
	
	void mate(uint16_t numMatingCouples);
	
	/**
	 * Transform the given fitness map to start at 0 and only have positive values
	 */
	void transformCurrentGenFitnessMap();

//	void cleanCurrentGeneration();

	const Control & ctrl;
	const Evaluator * const evaluator;
};

typedef std::multiset<Chromosome, Population::ChromosomeComparator> SortedChromosomes;

#endif
