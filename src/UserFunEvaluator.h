//
//  UserFunEvaluator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#ifndef GenAlgPLS_UserFunEvaluator_h
#define GenAlgPLS_UserFunEvaluator_h

#include "config.h"

#include <RcppArmadillo.h>
#include "Evaluator.h"
#include "Chromosome.h"

class UserFunEvaluator : public Evaluator {
public:
	UserFunEvaluator(Rcpp::Function const &userFun, const VerbosityLevel &verbosity) : Evaluator(verbosity), userFun(userFun) {};
//	~UserFunEvaluator();

	double evaluate(Chromosome &ch) const;

private:
	const Rcpp::Function userFun;
};

#endif
