//
//  Control.h
//
//

#ifndef GenAlgPLS_Control_h
#define GenAlgPLS_Control_h

#include "config.h"

#include <iostream>
#include <vector>
#include <RcppArmadillo.h>

enum VerbosityLevel {
	OFF = 0,
	ON,
	MORE_VERBOSE,
	DEBUG_VERBOSE,
	FULLY_VERBOSE
};

enum CrossoverType {
	SINGLE = 0,
	RANDOM = 1
};

class Control {
public:
	Control(const uint16_t chromosomeSize,
			const uint16_t popSize,
			const uint16_t numGenerations,
			const uint16_t elitism,
			const uint16_t minVariables,
			const uint16_t maxVariables,
			const uint16_t maxMatingTries,
			const double mutationProbability,
			const uint16_t numThreads,
			const uint16_t maxDuplicateEliminationTries,
			const enum CrossoverType crossover,
			const enum VerbosityLevel verbosity) : chromosomeSize(chromosomeSize), populationSize(popSize), numGenerations(numGenerations), elitism(elitism), minVariables(minVariables), maxVariables(maxVariables),	maxMatingTries(maxMatingTries), mutationProbability(mutationProbability), numThreads(numThreads), maxDuplicateEliminationTries(maxDuplicateEliminationTries), crossover(crossover), verbosity(verbosity) {};

	const uint16_t chromosomeSize;
	const uint16_t populationSize;
	const uint16_t numGenerations;
	const uint16_t elitism;
	const uint16_t minVariables;
	const uint16_t maxVariables;
	const uint16_t maxMatingTries;
	const double mutationProbability;
	const uint16_t numThreads;
	const uint16_t maxDuplicateEliminationTries;
	const enum CrossoverType crossover;
	const enum VerbosityLevel verbosity;

	friend std::ostream& operator<<(std::ostream &os, const Control &ctrl) {
		os << "Chromosome size: " << ctrl.chromosomeSize << std::endl
		<< "Population size: " << ctrl.populationSize << std::endl
		<< "Number of generations: " << ctrl.numGenerations << std::endl
		<< "Number of elite chromosomes to keep: " << ctrl.elitism << std::endl
		<< "Number of variables set: " << ctrl.minVariables << " to " << ctrl.maxVariables << std::endl
		<< "Maximum number of tries for mating: " << ctrl.maxMatingTries << std::endl
		<< "Mutation probability: " << ctrl.mutationProbability << std::endl
		<< "Maximum number of tries to eliminate duplicates: " << ctrl.maxDuplicateEliminationTries << std::endl
		<< "Crossover-type: " << ((ctrl.crossover == SINGLE) ? "Single" : "Random") << std::endl
		<< "Number of threads: " << ctrl.numThreads << std::endl
		<< "Verbosity Level: " << ctrl.verbosity << std::endl;
		
		return os;
		};
};

#endif
