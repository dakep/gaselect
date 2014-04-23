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
	PLSEvaluator(PLS* pls, uint16_t numReplications, uint16_t maxNComp, const std::vector<uint32_t> &seed, VerbosityLevel verbosity,
	uint16_t innerSegments, uint16_t outerSegments = 1, double testSetSize = 0.0);

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
	const uint16_t outerSegments;
	const uint16_t innerSegments;
	const arma::uword nrows;
	const bool cloned;

	PLS *pls;
	uint16_t maxNComp;
	std::vector<arma::uvec> segmentation;

	PLSEvaluator(const PLSEvaluator &other);

	/**
	 * Estimate the SEP
	 */
	double estSEP(uint16_t maxNComp);

	void initSegmentation(double testSetSize, const std::vector<uint32_t> &seed);

	/**
	 * Helper class for online calculation of the standard deviation (and optionally the sum)
	 */
	class OnlineStddev {
	public:
		OnlineStddev(uint16_t dim = 1) : dim(dim),
			meanVec(dim, 0.0), M2(dim, 0.0), counter(0) {}

		void reset() {
			this->counter = 0;
			std::memset(&(this->meanVec[0]), 0.0, this->dim * sizeof(this->meanVec[0]));
			std::memset(&(this->M2[0]), 0.0, this->dim * sizeof(this->M2[0]));
		};

		void update(arma::vec samples, uint16_t dim = 0) {
			for(uint16_t i = 0; i < samples.n_elem; ++i) {
				this->update(samples[i], dim);
			}
		};

		void update(double sample, uint16_t dim = 0) {
			double delta = sample - this->meanVec[dim];
			this->meanVec[dim] = this->meanVec[dim] + delta / (++this->counter);
			this->M2[dim] = this->M2[dim] + delta * (sample - this->meanVec[dim]);
		};

		double stddev(uint16_t dim = 0) {
			return std::sqrt(this->M2[dim] / (this->counter - 1));
		};

		double mean(uint16_t dim = 0) {
			return this->meanVec[dim];
		};

	private:
		const uint16_t dim;
		std::vector<double> meanVec;
		std::vector<double> M2;
		uint16_t counter;
	};
};


#endif
