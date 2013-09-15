//
//  GenAlg.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 05.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef GenAlgPLS_GenAlg_h
#define GenAlgPLS_GenAlg_h

#include "config.h"

#include <RcppArmadillo.h>

#ifdef TIMING_BENCHMARK

#include <sys/time.h>
#include <sys/types.h>

#endif

enum EvaluatorClass {
	USER = 0,
	PLS_EVAL = 1,
	LM = 2
};

/**
 * arguments:
 *	control ... A R list with following entries:
 *		uint16_t chromosomeSize ... The size of the chromosome (most times equal to the number of columns of X) (> 0)
 *		uint16_t populationSize ... Number of indivudual chromosomes in the population (i.e. per generation) (> 0)
 *		uint16_t numGenerations ... The number of generations to generate (> 0)
 *		uint16_t minVariables ... The minimum number of variables in a subset
 *		uint16_t maxVariables ... The maximum number of variables in a subset
 *		uint16_t maxMatingTries ... The maximum number of tries to get better children than parents
 *		uint16_t elitism ... The number of "elite" chromosomes to keep accross all generations (>= 0)
 *		double mutationProb ... The probability of using a new variable (0 <= onesRatio < 1)
 *		uint16_t numThreads ... The maximum number of threads to spawn
 *		VerbosityLevel verbosity ... Level of verbosity
 *		EvaluatorClass evaluatorClass ... The evaluator to use
 *		Rcpp::Function userEvalFunction ... The function to be called for evaluating the fitness of a chromosome
 *		PLSMethod plsMethod ... PLS method to use in internal evaluation
 *		uint16_t numReplications ... Number of replications in the internal evaluation procedure (the variable subset is evaluted with CV numReplication times and the mean fitness is returned) (> 0)
 *		uint16_t numSegments ... Number of CV segments used in the internal evaluation method (> 0)
 *		int statistic ... The statistic the LM Evaluator should use
 *	X ... A numeric matrix with dimensions n x p (optional - only needed if using internal evaluation methods)
 *	y ... A numeric vector with length n (optional - only needed if using internal evaluation methods)
 *	seed ... An integer (uint32_t) with the initial seed
 */
RcppExport SEXP genAlgPLS(SEXP control, SEXP X, SEXP y, SEXP seed);

/**
 * evaluate the given data with the given evaluator
 * arguments:
 *	evaluator ... A R list with following entries
 *		EvaluatorClass evaluatorClass ... The evaluator to use
 *		Rcpp::Function userEvalFunction ... The function to be called for evaluating the fitness of a chromosome
 *		PLSMethod plsMethod ... PLS method to use in internal evaluation
 *		uint16_t numReplications ... Number of replications in the internal evaluation procedure (the variable subset is evaluted with CV numReplication times and the mean fitness is returned) (> 0)
 *		uint16_t numSegments ... Number of CV segments used in the internal evaluation method (> 0)
 *		int statistic ... The statistic the LM Evaluator should use
 *
 *	X ... A numeric matrix with dimensions n x p
 *	y ... A numeric vector with length n
 *	seed ... An integer (uint32_t) with the initial seed
 */
RcppExport SEXP evaluate(SEXP evaluator, SEXP X, SEXP y, SEXP seed);

// RcppExport SEXP simpls(SEXP X, SEXP Y, SEXP ncomp, SEXP newX);
//
// RcppExport SEXP evalTest(SEXP X, SEXP Y, SEXP numReplications, SEXP numSegments);

#endif
