//
//  PLSEvaluator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#ifndef GenAlgPLS_PLSEvaluator_h
#define GenAlgPLS_PLSEvaluator_h

#include "config.h"

#include <stdexcept>
#include <vector>
#include <algorithm>
#include <RcppArmadillo.h>

#include "RNG.h"
#include "Evaluator.h"
#include "Chromosome.h"
#include "PLS.h"

class PLSEvaluator : public Evaluator {
public:
	PLSEvaluator(PLS* pls, const uint16_t numReplications, const uint16_t numSegments, const std::vector<uint32_t> &seed, const VerbosityLevel verbosity);

	~PLSEvaluator() {
		if(this->cloned == true) {
			delete this->pls;
		}
	}
	
	double evaluate(Chromosome &ch) {
		arma::uvec columnSubset = ch.toColumnSubset();
		double fitness = this->evaluate(columnSubset);
		ch.setFitness(fitness);
		return fitness;
	};

	double evaluate(arma::uvec &columnSubset);
	
	Evaluator* clone() const;

#ifdef ENABLE_DEBUG_VERBOSITY
	static uint32_t counter;
#endif

private:
	const uint16_t numReplications;
	const uint16_t numSegments;
	const arma::uword nrows;
	const arma::uword segmentLength; // The length of the incomplete segments
	const uint16_t completeSegments; // The number of segments with `segmentLength` + 1 elements. If 0, all segments have `segmentLength` elements
	const std::vector<uint32_t> &seed;

	PLS *pls;
	bool cloned;

	/**
	 * Estimate the SEP
	 */
	double estSEP(uint16_t ncomp, std::vector<arma::uword> &rowNumbers);

	std::vector< std::vector<arma::uword> > shuffledRowNumbers;
	void initRowNumbers();
};


#endif
