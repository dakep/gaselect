//
//  TruncatedGeomGenerator.h
//  GenAlgTest
//
//  Created by David Kepplinger on 18.05.2013.
//
//

#ifndef GenAlgPLS_TruncatedGeomGenerator_h
#define GenAlgPLS_TruncatedGeomGenerator_h

#include <RcppArmadillo.h>
#include <inttypes.h>

#include "UnifGenerator__0__1.h"

class TruncatedGeomGenerator : public Rcpp::Generator<uint16_t> {
public:
	TruncatedGeomGenerator(const double p) : prob(p), commonDenominator(log1p(-p)) {
	}

	/*
	 * Uses the inversion method for generating a truncated exponential variate and
	 * than returning the rounded-down value which is truncated geometrically distributed
	 * cutoff must be > 0 otherwise !
	 */
	inline uint16_t operator()(const uint16_t cutoff) const {
		return (uint16_t) (log1p(- this->unifGen() * (1. - R_pow_di(1. - this->prob, cutoff))) / this->commonDenominator);
	}

private:
	const Rcpp::stats::UnifGenerator__0__1 unifGen;
	const double prob;
	const double commonDenominator;
};

#endif
