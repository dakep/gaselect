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
	Population(const Control &ctrl, ::Evaluator &evaluator, SynchronizedUnifGenerator__0__1& unifGen);
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

#ifdef HAVE_PTHREAD_H
	class ThreadingError : public std::runtime_error {
	public:
		ThreadingError(const char* what) : std::runtime_error(what) {};
		virtual ~ThreadingError() throw() {};
	};
#endif
	
private:	
	std::multiset<Chromosome, Population::ChromosomeComparator> elite;
	std::vector<Chromosome> currentGeneration;
	std::vector<Chromosome> nextGeneration;
	std::vector<double> currentGenFitnessMap;
	std::vector<double> nextGenFitnessMap;
	double sumCurrentGenFitness;
	double minCurrentGenFitness; // Minimum fitness value in the current generation
	double minEliteFitness;

	/**
	 * Pick a chromosome from the current generation at random
	 * where the probability to pick a chromosome is taken from
	 * the currentGenFitnessMap
	 */
	Chromosome &drawChromosomeFromCurrentGeneration(SynchronizedUnifGenerator__0__1& unifGen);
	void addChromosomeToElite(Chromosome &ch);

	std::ostream& printChromosomeFitness(std::ostream &os, Chromosome &ch);
	
	void mate(uint16_t numMatingCouples, ::Evaluator& evaluator, SynchronizedUnifGenerator__0__1& unifGen, bool checkUserInterrupt = true);
	
	/**
	 * Transform the given fitness map to start at 0 and only have positive values
	 */
	void transformCurrentGenFitnessMap();

//	void cleanCurrentGeneration();

	const Control& ctrl;
	::Evaluator& evaluator;
	SynchronizedUnifGenerator__0__1& unifGen;

#ifdef HAVE_PTHREAD_H	
	struct ThreadArgsWrapper {
		Population* popObj;
		Evaluator* evalObj;
		SynchronizedUnifGenerator__0__1* unifGen;
		uint16_t numMatingCouples;
	};
	
	static void* matingThreadStart(void* obj);
	
	void runMating(uint16_t numMatingCoupls, ::Evaluator& evaluator, SynchronizedUnifGenerator__0__1& unifGen);
	pthread_t* threads;
	
	/*
	 * Mutex and condition variables
	 */
	
	pthread_mutex_t printMutex;
	pthread_mutex_t updateMutex;
	pthread_mutex_t syncMutex;
	pthread_cond_t startMatingCond;
	pthread_cond_t allThreadsFinishedMatingCond;
	
	bool startMating;
	bool killThreads;
	bool allThreadsFinishedMating;
	
	uint16_t actuallySpawnedThreads;
	uint16_t numThreadsFinishedMating;
	
	void finishedMating();
	void waitForAllThreadsToFinishMating();
	
	void cancelAllThreads();
	
	void cleanupPthread();
#endif
};

typedef std::multiset<Chromosome, Population::ChromosomeComparator> SortedChromosomes;

#endif
