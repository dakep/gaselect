//
//  ThreadedPopulation.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//

#ifndef GenAlgPLS_ThreadedPopulation_h
#define GenAlgPLS_ThreadedPopulation_h
#include "config.h"

#ifdef HAVE_PTHREAD_H
#include <iostream>
#include <vector>
#include <set>

#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"
#include "Population.h"

#include "RNG.h"

class MultiThreadedPopulation : public Population {
	
public:
	MultiThreadedPopulation(const Control &ctrl, ::Evaluator &evaluator, RNG& rng);
	~MultiThreadedPopulation();
	
	/**
	 * Returns the last generation and the elite (if any)
	 * Invalid once the Population object is destroyed!
	 */
	
	class ThreadingError : public std::runtime_error {
	public:
		ThreadingError(const char* what) : std::runtime_error(what) {};
		virtual ~ThreadingError() throw() {};
	};
	
	void run();
private:
	struct ThreadArgsWrapper {
		MultiThreadedPopulation* popObj;
		Evaluator* evalObj;
		unsigned int seed;
		uint16_t numMatingCouples;
		uint16_t offset;
	};
	
	std::vector<Chromosome*> nextGeneration;
	std::vector<double> nextGenFitnessMap;
	double sumCurrentGenFitness;
	double minCurrentGenFitness; // Minimum fitness value in the current generation
	
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
		
	void mate(uint16_t numMatingCouples, ::Evaluator& evaluator, RNG& rng, uint16_t offset, bool checkUserInterrupt = true);
	
	/**
	 * Transform the given fitness map to start at 0 and only have positive values
	 */
	void transformCurrentGenFitnessMap();
	
	static void* matingThreadStart(void* obj);
	
	void runMating(uint16_t numMatingCoupls, ::Evaluator& evaluator, RNG& rng, uint16_t offset);
	void finishedMating();
	void waitForAllThreadsToFinishMating();
};


#endif
#endif
