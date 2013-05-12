//
//  GenAlg.cpp
//  GenAlgTest
//

#include "config.h"

#include <RcppArmadillo.h>
#include <set>
#include <inttypes.h>

#include "Chromosome.h"
#include "Population.h"
#include "Control.h"
#include "UserFunEvaluator.h"
#include "PLSEvaluator.h"


#include "GenAlg.h"

using namespace Rcpp;

SEXP genAlg(SEXP Scontrol, SEXP SX, SEXP Sy) {
	::Evaluator *eval;
	PLS *pls;
	uint8_t toFree = 0; // first bit is set ==> free eval; 2nd bit set ==> free pls
	bool useUserFunction = false;
BEGIN_RCPP

	Control ctrl;
	List control = List(Scontrol);

	// All checks are disabled and must be performed in the R code calling this script
	// Otherwise unexpected behaviour

	// All list elements accessed below HAVE to be set in the control list
	ctrl.setChromosomeSize(as<uint16_t>(control["chromosomeSize"]));
	ctrl.setPopulationSize(as<uint16_t>(control["populationSize"]));
	ctrl.setNumberOfGenerations(as<uint16_t>(control["numGenerations"]));
	ctrl.setElitism(as<uint16_t>(control["elitism"]));
	ctrl.setMutate0To1Probability(as<double>(control["mutationProb"]));
	ctrl.setOnesRatio(as<double>(control["onesRatio"]));
	ctrl.setVerbosity((VerbosityLevel) as<int>(control["verbosity"]));

	useUserFunction = as<bool>(control["useUserSuppliedFunction"]);
	if(useUserFunction) {
		eval = new UserFunEvaluator(as<Rcpp::Function>(control["userEvalFunction"]), ctrl.getVerbosity());
	} else {
		Rcpp::NumericMatrix XMat(SX);
		Rcpp::NumericMatrix YMat(Sy);
		arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false, true);
		arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false, true);
		PLSMethod method = (PLSMethod) as<int>(control["plsMethod"]);

		pls = PLS::getInstance(method, X, Y, false);
		toFree |= 2; // pls has to be freed
		
		eval = new PLSEvaluator(*pls, as<uint16_t>(control["numReplications"]), as<uint16_t>(control["numSegments"]), ctrl.getVerbosity());
	}
	toFree |= 1; // eval has to be freed

	if(ctrl.getVerbosity() == MORE_VERBOSE) {
		Rcout << ctrl << std::endl;		
	}
	
	Population pop(ctrl, *eval);
	pop.run();

	SortedChromosomes result = pop.getResult();

	Rcpp::LogicalMatrix retMatrix(ctrl.getChromosomeSize(), (const int) result.size());
	Rcpp::NumericVector retFitnesses((const int) result.size());
	uint16_t i = (uint16_t) result.size() - 1;
	
	for(SortedChromosomes::iterator it = result.begin(); it != result.end(); ++it, --i) {
		retFitnesses[i] = it->getFitness();
		retMatrix.column(i) = it->toLogicalVector();
	}

	delete eval; // must definitely be freed
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
	
	return R_NilValue;
}

SEXP simpls(SEXP Xs, SEXP Ys, SEXP numOfComp, SEXP newXs) {
BEGIN_RCPP	
	Rcpp::NumericMatrix XMat(Xs);
	Rcpp::NumericMatrix YMat(Ys);
	Rcpp::NumericMatrix newXMat(newXs);
	uint16_t ncomp = Rcpp::as<uint16_t>(numOfComp);
	
	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
	arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
	arma::mat newX(newXMat.begin(), newXMat.nrow(), newXMat.ncol(), false);

	PLSSimpls simpls(X, Y, true);
	simpls.fit(ncomp);

	return Rcpp::List::create(Rcpp::Named("coefficients") = simpls.getCoefficients(),
							  Rcpp::Named("fitted.values") = simpls.getFittedValues(),
							  Rcpp::Named("predicted") = simpls.predict(newX));
END_RCPP
}


RcppExport SEXP evalTest(SEXP Xs, SEXP Ys, SEXP numReplications, SEXP numSegments) {
BEGIN_RCPP
	Rcpp::NumericMatrix XMat(Xs);
	Rcpp::NumericMatrix YMat(Ys);
	arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
	arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
	PLSSimpls pls(X, Y, false);
	PLSEvaluator eval(pls, Rcpp::as<uint16_t>(numReplications), Rcpp::as<uint16_t>(numSegments), DEBUG_VERBOSE);

	arma::uvec colSubset(X.n_cols);
	
	for(uint16_t i = 0; i < X.n_cols; ++i) {
		colSubset[i] = i;
	}
	
	return Rcpp::wrap(eval.evaluate(colSubset));
END_RCPP
}

