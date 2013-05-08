//
//  UserFunEvaluator.cpp
//  GenAlgTest
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#include "UserFunEvaluator.h"

double UserFunEvaluator::evaluate(Chromosome &ch) const {
	SEXP rawFitness = this->userFun(Rcpp::wrap(ch.toLogicalVector()));
	if(!Rf_isNumeric(rawFitness)) {
		throw Rcpp::exception("Evaluation function has to return a numeric value", __FILE__, __LINE__);
	}
	
	double fitness = Rcpp::as<double>(rawFitness);
	ch.setFitness(fitness);
	return fitness;
}