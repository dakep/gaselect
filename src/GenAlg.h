//
//  GenAlg.h
//  GenAlgTest
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

/**
 * arguments:
 *	control ... A R list with following entries:
 *		uint16_t chromosomeSize ... The size of the chromosome (most times equal to the number of columns of X) (> 0)
 *		uint16_t populationSize ... Number of indivudual chromosomes in the population (i.e. per generation) (> 0)
 *		uint16_t numGenerations ... The number of generations to generate (> 0)
 *		uint16_t elitism ... The number of "elite" chromosomes to keep accross all generations (>= 0)
 *		double mutationProb ... The probability of using a new variable (0 <= onesRatio < 1)
 *		double onesRatio ... The approximate ratio of 1's in the chromosome (is actually distributed binomial with parameter onesRatio) (0 < onesRatio < 1)
 *		VerbosityLevel verbosity ... Level of verbosity
 *
 *		bool useUserSuppliedFunction ... If true, a user specified function is used to evaluate the fitness of a chromosome
 *		Rcpp::Function userEvalFunction ... The function to be called for evaluating the fitness of a chromosome
 *		PLSMethod plsMethod ... PLS method to use in internal evaluation
 *		uint16_t numReplications ... Number of replications in the internal evaluation procedure (the variable subset is evaluted with CV numReplication times and the mean fitness is returned) (> 0)
 *		uint16_t numSegments ... Number of CV segments used in the internal evaluation method (> 0)
 *	X ... A numeric matrix with dimensions n x p (optional - only needed if using internal evaluation methods)
 *	y ... A numeric vector with length n (optional - only needed if using internal evaluation methods)
 */
RcppExport SEXP genAlgPLS(SEXP control, SEXP X, SEXP y);

// RcppExport SEXP simpls(SEXP X, SEXP Y, SEXP ncomp, SEXP newX);
//
// RcppExport SEXP evalTest(SEXP X, SEXP Y, SEXP numReplications, SEXP numSegments);

#endif
