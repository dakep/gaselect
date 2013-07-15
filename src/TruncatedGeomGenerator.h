//
//  TruncatedGeomGenerator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 18.05.2013.
//
//

#ifndef GenAlgPLS_TruncatedGeomGenerator_h
#define GenAlgPLS_TruncatedGeomGenerator_h

#include "config.h"

#include <RcppArmadillo.h>

#include "SynchronizedUnifGenerator__0__1.h"

#if RCPP_VERSION < Rcpp_Version(0, 10, 1)
class TruncatedGeomGenerator : public Rcpp::Generator<false, uint16_t> {
#else
class TruncatedGeomGenerator : public Rcpp::Generator<uint16_t> {
#endif
public:
	TruncatedGeomGenerator(const double p) : unifGen(), prob(p), commonDenominator(log1p(-p)) {
	}

	/*
	 * Uses the inversion method for generating a truncated exponential variate and
	 * than returning the rounded-down value which is truncated geometrically distributed
	 * cutoff must be > 0 otherwise undefined behaviour!
	 */
	inline uint16_t operator()(const uint16_t cutoff) const {
		return (uint16_t) (log1p(- TruncatedGeomGenerator::unifGen() * (1. - R_pow_di(1. - this->prob, cutoff + 1))) / this->commonDenominator);
	}

private:
	SynchronizedUnifGenerator__0__1 unifGen;

	const double prob;
	const double commonDenominator;
};

#endif
