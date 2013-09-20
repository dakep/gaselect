//
//  GenAlg.cpp
//  GenAlgPLS
//

#include "config.h"

#include <RcppArmadillo.h>
#include <set>

#include "Chromosome.h"
#include "Control.h"
#include "PLS.h"
#include "UserFunEvaluator.h"
#include "PLSEvaluator.h"
#include "LMEvaluator.h"
#include "SingleThreadPopulation.h"
#include "RNG.h"

#ifdef HAVE_PTHREAD_H
#include "MultiThreadedPopulation.h"
#endif

#include "GenAlg.h"

using namespace Rcpp;

SEXP genAlgPLS(SEXP Scontrol, SEXP SX, SEXP Sy, SEXP Sseed) {
	::Evaluator *eval;
	PLS *pls;
	Population *pop;
	uint8_t toFree = 0; // first bit is set ==> free eval; 2nd bit set ==> free pls; 3rd bit set ==> free globalUnifGen; 4th bit set ==> free pop
BEGIN_RCPP
	List control = List(Scontrol);
	uint32_t singleSeed = as<uint32_t>(Sseed);
	std::vector<uint32_t> seed;
	uint16_t numThreads = as<uint16_t>(control["numThreads"]);
	VerbosityLevel verbosity = (VerbosityLevel) as<int>(control["verbosity"]);
	EvaluatorClass evalClass = (EvaluatorClass) as<int>(control["evaluatorClass"]);
	
	if(numThreads > 1) {
#ifdef HAVE_PTHREAD_H
		if(evalClass == USER) {
			Rcout << "Warning: Multithreading is not available when using a user supplied function for evaluation" << std::endl;
		}
		
		if(verbosity >= MORE_VERBOSE) {
			verbosity = MORE_VERBOSE;
		}
#else
		Rcout << "Warning: Threads are not supported on this system" << std::endl;
		numThreads = 1;
#endif
	} else if(numThreads < 1) {
		numThreads = 1;
	}
	
	// All checks are disabled and must be performed in the R code calling this script
	// Otherwise unexpected behaviour
	Control ctrl(as<uint16_t>(control["chromosomeSize"]),
				 as<uint16_t>(control["populationSize"]),
				 as<uint16_t>(control["numGenerations"]),
				 as<uint16_t>(control["elitism"]),
				 as<uint16_t>(control["minVariables"]),
				 as<uint16_t>(control["maxVariables"]),
				 as<uint16_t>(control["maxMatingTries"]),
				 as<double>(control["mutationProb"]),
				 numThreads,
				 as<uint8_t>(control["cutoffQuantile"]),
				 (CrossoverType) as<int>(control["crossover"]),
				 verbosity);

	/*
	 * Generate a common seed for the Population and the PLSEvaluator objects
	 */
	RNG rng(singleSeed);
	seed.reserve(RNG::SEED_SIZE);
	for(uint32_t i = 0; i < RNG::SEED_SIZE; ++i) {
		seed.push_back(rng());
	}
	
	switch(evalClass) {
		case USER: {
			eval = new UserFunEvaluator(as<Rcpp::Function>(control["userEvalFunction"]), ctrl.verbosity);
			break;
		}
		case PLS_EVAL: {
			Rcpp::NumericMatrix XMat(SX);
			Rcpp::NumericMatrix YMat(Sy);
			arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
			arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
			PLSMethod method = (PLSMethod) as<int>(control["plsMethod"]);
			
			pls = PLS::getInstance(method, X, Y, false);
			toFree |= 2; // pls has to be freed
			
			eval = new PLSEvaluator(pls, as<uint16_t>(control["numReplications"]), as<uint16_t>(control["numSegments"]), seed, ctrl.verbosity);

			break;
		}
		case LM: {
			Rcpp::NumericMatrix XMat(SX);
			Rcpp::NumericMatrix YMat(Sy);
			arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
			arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
			arma::colvec y = Y.col(0);
			
			LMEvaluator::Statistic stat = (LMEvaluator::Statistic) as<int>(control["statistic"]);
			eval = new LMEvaluator(X, y, stat, verbosity);
			
			break;
		}
		default:
			break;
	}
	
	toFree |= 1; // eval has to be freed

	if(ctrl.verbosity >= MORE_VERBOSE) {
		Rcout << ctrl << std::endl;
	}

#ifdef HAVE_PTHREAD_H
	try {
		if(numThreads > 1) {
			pop = new MultiThreadedPopulation(ctrl, *eval, seed);
		} else {
			pop = new SingleThreadPopulation(ctrl, *eval, seed);
		}
		toFree |= 8;
		pop->run();
	} catch(MultiThreadedPopulation::ThreadingError& te) {
		if(ctrl.verbosity >= DEBUG_VERBOSE) {
			throw te;
		} else {
			throw Rcpp::exception("Multithreading could not be initialized. Set numThreads to 0 to avoid this problem.", __FILE__, __LINE__);
		}
	}
#else
	pop = new SingleThreadPopulation(ctrl, *eval);
	toFree |= 8;
	pop->run();
#endif
	SortedChromosomes result = pop->getResult();

	Rcpp::LogicalMatrix retMatrix(ctrl.chromosomeSize, (const int) result.size());
	Rcpp::NumericVector retFitnesses((const int) result.size());
	uint16_t i = (uint16_t) result.size() - 1;

	for(SortedChromosomes::iterator it = result.begin(); it != result.end(); ++it, --i) {
		retFitnesses[i] = it->getFitness();
		retMatrix.column(i) = it->toLogicalVector();
	}

	delete eval;
	
	if((toFree & 8) > 0) {
		delete pop;
	}

	if((toFree & 2) > 0) {
		delete pls;
	}

	return Rcpp::List::create(Rcpp::Named("subsets") = retMatrix,
							  Rcpp::Named("fitness") = retFitnesses);
VOID_END_RCPP
	if((toFree & 1) > 0) {
		delete eval;
	}
	if((toFree & 2) > 0) {
		delete pls;
	}
	if((toFree & 8) > 0) {
		delete pop;
	}

	return R_NilValue;
}



/**
 *
 */
SEXP evaluate(SEXP Sevaluator, SEXP SX, SEXP Sy, SEXP Sseed) {
	double fitness = 0.0;
	::Evaluator *eval;
	PLS *pls;
	uint8_t toFree = 0; // first bit is set ==> free eval; 2nd bit set ==> free pls; 3rd bit set ==> free globalUnifGen; 4th bit set ==> free pop
	
	BEGIN_RCPP
	List evaluator = List(Sevaluator);
	Rcpp::NumericMatrix XMat(SX);
	Rcpp::NumericMatrix YMat(Sy);
	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
	arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
	EvaluatorClass evalClass = (EvaluatorClass) as<int>(evaluator["evaluatorClass"]);
	std::vector<uint32_t> seed;

	
	switch(evalClass) {
		case USER: {
			eval = new UserFunEvaluator(as<Rcpp::Function>(evaluator["userEvalFunction"]), OFF);
			break;
		}
		case PLS_EVAL: {
			PLSMethod method = (PLSMethod) as<int>(evaluator["plsMethod"]);
			RNG rng(as<uint32_t>(Sseed));

			seed.reserve(RNG::SEED_SIZE);
			for (uint32_t i = 0; i < RNG::SEED_SIZE; ++i) {
				seed.push_back(rng());
			}
			
			pls = PLS::getInstance(method, X, Y, false);
			toFree |= 2; // pls has to be freed
		
			eval = new PLSEvaluator(pls, as<uint16_t>(evaluator["numReplications"]), as<uint16_t>(evaluator["numSegments"]), seed, OFF);
			
			break;
		}
		case LM: {
			arma::colvec y = Y.col(0);
			
			LMEvaluator::Statistic stat = (LMEvaluator::Statistic) as<int>(evaluator["statistic"]);
			eval = new LMEvaluator(X, y, stat, OFF);
			
			break;
		}
		default:
			break;
	}
	
	toFree |= 1; // eval has to be freed
	
	arma::uvec allCols(X.n_cols);

	for(arma::uword i = 0; i < allCols.n_elem; ++i) {
		allCols[i] = i;
	}
	
	fitness = eval->evaluate(allCols);
	
	if((toFree & 2) > 0) {
		delete pls;
	}
	
	return Rcpp::wrap(fitness);
	
	VOID_END_RCPP
	if((toFree & 1) > 0) {
		delete eval;
	}
	if((toFree & 2) > 0) {
		delete pls;
	}
	
	return R_NilValue;
}

//RcppExport SEXP WELL19937a(SEXP Sn, SEXP Smin, SEXP Smax, SEXP SnStreams, SEXP Sseed) {
//	uint32_t n = as<uint32_t>(Sn);
//	uint32_t nStreams = as<uint32_t>(SnStreams);
//	uint32_t seed = as<uint32_t>(Sseed);
//	double min = as<double>(Smin);
//	double max = as<double>(Smax);
//	uint16_t row = 0;
//	
//	RNG rng(seed);
//	
//	Rcpp::NumericMatrix retMat(n, nStreams);
//
//	for(uint16_t col = 0; col < nStreams; ++col) {
//		for(row = 0; row < n; ++row) {
//			retMat.column(col)[row] = rng(min, max);
//		}
//	}
//
//	Rcpp::Rcout << retMat.rows() << " x " << retMat.cols() << std::endl;
//	
//	return Rcpp::wrap(retMat);
//}

// SEXP simpls(SEXP Xs, SEXP Ys, SEXP numOfComp, SEXP newXs) {
// BEGIN_RCPP
// 	Rcpp::NumericMatrix XMat(Xs);
// 	Rcpp::NumericMatrix YMat(Ys);
// 	Rcpp::NumericMatrix newXMat(newXs);
// 	uint16_t ncomp = Rcpp::as<uint16_t>(numOfComp);
//
// 	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
// 	arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
// 	arma::mat newX(newXMat.begin(), newXMat.nrow(), newXMat.ncol(), false);
//
// 	PLSSimpls simpls(X, Y, true);
// 	simpls.fit(ncomp);
//
// 	return Rcpp::List::create(Rcpp::Named("coefficients") = simpls.getCoefficients(),
// 							  Rcpp::Named("fitted.values") = simpls.getFittedValues(),
// 							  Rcpp::Named("predicted") = simpls.predict(newX));
// END_RCPP
// }
//
//
// RcppExport SEXP evalTest(SEXP Xs, SEXP Ys, SEXP numReplications, SEXP numSegments) {
// BEGIN_RCPP
// 	Rcpp::NumericMatrix XMat(Xs);
// 	Rcpp::NumericMatrix YMat(Ys);
// 	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
// 	arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
// 	PLSSimpls pls(X, Y, false);
// 	PLSEvaluator eval(pls, Rcpp::as<uint16_t>(numReplications), Rcpp::as<uint16_t>(numSegments), DEBUG_VERBOSE);
//
// 	arma::uvec colSubset(X.n_cols);
//
// 	for(uint16_t i = 0; i < X.n_cols; ++i) {
// 		colSubset[i] = i;
// 	}
//
// 	return Rcpp::wrap(eval.evaluate(colSubset));
// END_RCPP
// }
//
