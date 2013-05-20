//
//  Control.h
//  
//

#ifndef GenAlgPLS_Control_h
#define GenAlgPLS_Control_h

#include "config.h"

#include <iostream>
#include <vector>
#include <inttypes.h>
#include <RcppArmadillo.h>
#include <Rcpp/stats/random/runif.h>

enum VerbosityLevel {
	OFF = 0,
	ON,
	MORE_VERBOSE,
	DEBUG_VERBOSE
};

class Control {

public:
	Control(const uint16_t chromosomeSize, const uint16_t popSize, const uint16_t numGenerations, const uint16_t elitism, const uint16_t minVariables, const uint16_t maxVariables, const double mutationProbability, const enum VerbosityLevel verbosity);
	const uint16_t chromosomeSize;
	const uint16_t populationSize;
	const uint16_t numGenerations;
	const uint16_t elitism;
	const uint16_t minVariables;
	const uint16_t maxVariables;
	const double mutationProbability;
	const enum VerbosityLevel verbosity;
	
	/*
	 * The values in this vector may be shuffled but the size must not be changed!
	 */
	std::vector<uint16_t> variablePositionPopulation;
	
	friend std::ostream& operator<<(std::ostream &os, const Control &ctrl);
	
	/**
	 * Convenience method to shuffle the variable-position population
	 * optionally only the first `length` positions in the population
	 * are random (though all positions are considered)
	 */
	void shuffleVariablePositionPopulation(uint16_t length = 0);

private:
	const Rcpp::stats::UnifGenerator__0__1 unifGen;
};

#endif
