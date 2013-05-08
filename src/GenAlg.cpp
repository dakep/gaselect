//
//  GenAlg.cpp
//  GenAlgTest
//

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
//SEXP genAlg(SEXP chromosomeSize, SEXP populationSize, SEXP numGenerations, SEXP onesRatio, SEXP mutationProb, SEXP elitism, SEXP verbosity, SEXP evalFunction) {
	
BEGIN_RCPP	
	Control ctrl;
	List control = List(Scontrol);
	bool useUserFunction = as<bool>(control["useUserSuppliedFunction"]);
	::Evaluator *eval;
	PLS *pls;
	
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
	
	if(useUserFunction) {
		eval = new UserFunEvaluator(as<Rcpp::Function>(control["userEvalFunction"]), ctrl.getVerbosity());
	} else {
		Rcpp::NumericMatrix XMat(SX);
		Rcpp::NumericMatrix YMat(Sy);
		arma::mat X(XMat.begin(), XMat.nrow(), XMat.ncol(), false);
		arma::mat Y(YMat.begin(), YMat.nrow(), YMat.ncol(), false);
		PLSMethod method = (PLSMethod) as<int>(control["plsMethod"]);

		pls = PLS::getInstance(method, X, Y, false);
		
		eval = new PLSEvaluator(*pls, as<uint16_t>(control["numReplications"]), as<uint16_t>(control["numSegments"]), ctrl.getVerbosity());
	}
	
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
		retFitnesses[i] = (*it)->getFitness();
		retMatrix.column(i) = (*it)->toLogicalVector();
	}

	delete eval;
	if(!useUserFunction) {
		delete pls;
	}
	
	return Rcpp::List::create(Rcpp::Named("subsets") = retMatrix,
							  Rcpp::Named("fitness") = retFitnesses);
END_RCPP
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
	PLSEvaluator eval(pls, Rcpp::as<uint16_t>(numReplications), Rcpp::as<uint16_t>(numSegments), DEBUG);
	
	arma::uvec colSubset(X.n_cols);
	
	for(uint16_t i = 0; i < X.n_cols; ++i) {
		colSubset[i] = i;
	}
	
	return Rcpp::wrap(eval.evaluate(colSubset));
END_RCPP
}

