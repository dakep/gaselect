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
#include <stdexcept>
#include <iostream>
#include <streambuf>
#include <vector>
#include <set>
#include <string>
#include <utility>

#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"
#include "Population.h"

#include "RNG.h"

class MultiThreadedPopulation : public Population {
	
public:
	MultiThreadedPopulation(const Control &ctrl, ::Evaluator &evaluator, const std::vector<uint32_t> &seed);
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
	class SafeRstreambuf : public std::streambuf {
	public:
		SafeRstreambuf();
		~SafeRstreambuf();
		void realFlush();
	protected:
		virtual std::streamsize xsputn(const char* s, std::streamsize n);
		virtual int	overflow(int c = EOF);
		virtual int sync();
	private:
		std::string buffer;
		pthread_mutex_t printMutex;
	};
	
	class SafeRostream : public std::ostream {
	public:
		SafeRostream() : std::ostream(new SafeRstreambuf()) {
			this->buf = static_cast<SafeRstreambuf*>(this->rdbuf());
		};
		
		void realFlush() {
			if(this->buf != NULL) {
				this->buf->realFlush();
			}
		};

		~SafeRostream() {
			if(this->buf != NULL) {
				delete this->buf;
				this->buf = NULL;
			}
		}
   private:
	   SafeRstreambuf* buf;
   };
	
	struct ThreadArgsWrapper {
		MultiThreadedPopulation* popObj;
		Evaluator* evalObj;
		uint32_t seed;
		uint16_t numChildren;
		uint16_t offset;
	};
	
	ChromosomeVec nextGeneration;
	double sumCurrentGenFitness;
	double minCurrentGenFitness; // Minimum fitness value in the current generation
	
	/*
	 * Mutex and condition variables
	 */
	pthread_mutex_t syncMutex;
	pthread_cond_t startMatingCond;
	pthread_cond_t allThreadsFinishedMatingCond;
	
	bool startMating;
	bool killThreads;
	bool allThreadsFinishedMating;
	
	SafeRostream Rout;
	
	uint16_t actuallySpawnedThreads;
	uint16_t numThreadsFinishedMating;
		
	void mate(uint16_t numChildren, ::Evaluator& evaluator, RNG& rng, uint16_t offset, bool checkUserInterrupt = true);
	
	static void* matingThreadStart(void* obj);
	
	void runMating(uint16_t numMatingCoupls, ::Evaluator& evaluator, RNG& rng, uint16_t offset);
	void waitForAllThreadsToFinishMating();
	
	class OrderChromosomePtr : public std::binary_function<Chromosome*, Chromosome*, bool> {
	public:
		bool operator()(const Chromosome* const ch1, const Chromosome* const ch2) {
			return (ch1->getFitness() < ch2->getFitness());
		};
	};
};


#endif
#endif
