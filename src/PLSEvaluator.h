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

#include <RcppArmadillo.h>
#include "RNG.h"

#include "Evaluator.h"
#include "Chromosome.h"
#include "PLS.h"

class PLSEvaluator : public Evaluator {
public:
	PLSEvaluator(PLS* pls, const uint16_t numReplications, const uint16_t numSegments, const VerbosityLevel verbosity, RNG* rng) :
	Evaluator(verbosity), numReplications(numReplications), numSegments(numSegments),
	nrows(pls->getNumberOfObservations()), segmentLength(nrows / numSegments),
	completeSegments(nrows % numSegments), pls(pls), rng(rng), cloned(false)
	{
		if(pls->getNumberOfResponseVariables() > 1) {
			throw Rcpp::exception("PLS evaluator only available for models with 1 response variable", __FILE__, __LINE__);
		}
	};

	~PLSEvaluator() {
		if(this->cloned == true) {
			delete this->pls;
		}
	}
	
	double evaluate(Chromosome &ch) const {
		arma::uvec columnSubset = ch.toColumnSubset();
		double fitness = this->evaluate(columnSubset);
		ch.setFitness(fitness);
		return fitness;
	};

	double evaluate(arma::uvec &colSubset) const;

	void setRNG(RNG* rng) {
		this->rng = rng;
	};
	
	Evaluator* clone() const;

private:
	const uint16_t numReplications;
	const uint16_t numSegments;
	const arma::uword nrows;
	const arma::uword segmentLength; // The length of the incomplete segments
	const uint16_t completeSegments; // The number of segments with `segmentLength` + 1 elements. If 0, all segments have `segmentLength` elements

	PLS *pls;
	RNG *rng;

	bool cloned;

	/**
	 * Estimate the SEP
	 */
	double estSEP(uint16_t ncomp, arma::uvec &rowNumbers) const;

	arma::uvec initRowNumbers() const;
};


#endif
