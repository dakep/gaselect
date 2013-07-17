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
#include "UnifGenerator_0_1.h"

class TruncatedGeomGenerator {
public:
	TruncatedGeomGenerator(const double p) : prob(p), commonDenominator(log1p(-p)) {
	}

	/*
	 * Uses the inversion method for generating a truncated exponential variate and
	 * than returning the rounded-down value which is truncated geometrically distributed
	 * cutoff must be > 0 otherwise undefined behaviour!
	 */
	inline uint16_t operator()(const uint16_t cutoff, UnifGenerator_0_1& unifGen) const {
		return (uint16_t) (log1p(- unifGen() * (1. - R_pow_di(1. - this->prob, cutoff + 1))) / this->commonDenominator);
	}

private:
	const double prob;
	const double commonDenominator;
};

#endif
