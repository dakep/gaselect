//
//  PLSEvaluator.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#include "config.h"

#include <algorithm>
#include "Logger.h"
#include "PLSEvaluator.h"

#ifdef ENABLE_DEBUG_VERBOSITY
#define IF_DEBUG(expr) if(this->verbosity == DEBUG_EVAL || this->verbosity == DEBUG_ALL) { expr; }
#else
#define IF_DEBUG(expr)
#endif

#ifdef ENABLE_DEBUG_VERBOSITY
uint32_t PLSEvaluator::counter = 0;
#endif

PLSEvaluator::PLSEvaluator(PLS* pls, uint16_t numReplications, uint16_t maxNComp, const std::vector<uint32_t> &seed, VerbosityLevel verbosity,
	uint16_t innerSegments, uint16_t outerSegments, double testSetSize) :
	Evaluator(verbosity), numReplications(numReplications),
	outerSegments((outerSegments < 1) ? 1 : outerSegments),
	innerSegments((outerSegments <= 1 && testSetSize == 0.0) ? innerSegments - 1 : innerSegments),
	nrows(pls->getNumberOfObservations()), cloned(false), pls(pls), maxNComp(maxNComp)
{
	/* assert outerSegments > 0 */
	if(pls->getNumberOfResponseVariables() > 1) {
		throw std::invalid_argument("PLS evaluator only available for models with 1 response variable");
	}

	if(outerSegments > 1) {
		testSetSize = 1.0 / (double) outerSegments;
	}

	if(testSetSize < 0.0 || testSetSize >= 1.0) {
		throw std::invalid_argument("The test set size must be within the interval (0, 1)");
	}

	IF_DEBUG(GAout << "Test set size: " << testSetSize << " -- outerSegments: " << outerSegments << std::endl);

	this->initSegmentation(testSetSize, seed);
}

PLSEvaluator::PLSEvaluator(const PLSEvaluator &other) :
	Evaluator(other.verbosity), numReplications(other.numReplications), outerSegments(other.outerSegments),
	innerSegments(other.innerSegments), nrows(other.nrows), cloned(true), maxNComp(other.maxNComp), segmentation(other.segmentation)
{
	this->pls = other.pls->clone();
}

double PLSEvaluator::evaluate(arma::uvec &columnSubset) {
	if(columnSubset.n_elem == 0) {
		GAerr << GAerr.lock() << "Can not evaluate empty variable subset" << GAerr.unlock();
		throw Evaluator::EvaluatorException("Can not evaluate empty variable subset");
	}
#ifdef ENABLE_DEBUG_VERBOSITY
	++PLSEvaluator::counter;
#endif
	this->pls->viewSelectColumns(columnSubset);

	double sumSEP = this->estSEP(((this->maxNComp < columnSubset.n_elem) ? this->maxNComp : columnSubset.n_elem));

	IF_DEBUG(GAout << "EVALUATOR: Sum of SEP:" << std::endl << sumSEP << std::endl)
	return sumSEP;
}

/**
 * Initialize the row-segmentation for each replication and all segmentations
 */
inline void PLSEvaluator::initSegmentation(double testSetSize, const std::vector<uint32_t> &seed) {
	RNG rng(seed);
	ShuffledSet rowNumbers(this->nrows);
	arma::uword i, n = 0;

	if(testSetSize == 0.0 && this->outerSegments == 1) {
		testSetSize = 1.0 / (this->innerSegments + 1.0);
	}

	/*
	 * The size of the outer segment and the number of outer segments with one extra
	 * observation
	 */
	arma::uword outerSegmentLength = this->nrows * testSetSize;
	arma::uword outerSegmentLengthRem = this->nrows % outerSegmentLength;

	this->segmentation.reserve(2 * this->numReplications * (this->innerSegments + 1) * this->outerSegments);

	/*
	 * The size of the inner segment and the number of inner segments with one extra
	 * observation, if the smaller outer segment is used.
	 * If the larger outer segment is used, one inner segment with an extra observation
	 * won't get this extra observation. If no inner segment has an extra observation,
	 * one inner segment will have one observation less.
	 */
	arma::uword innerSegmentLength = (this->nrows - outerSegmentLength) / this->innerSegments;
	arma::uword innerSegmentLengthRem = (this->nrows - outerSegmentLength) % this->innerSegments;
	arma::uword innerSegmentLengthBigOuter = innerSegmentLength;

	if(outerSegmentLengthRem > 0 && innerSegmentLengthRem == 0) {
		--innerSegmentLengthBigOuter;
	}

	/*
	 * Update the maximal number of components, as innerSegmentLengthBigOuter is the minimal
	 * segment length.
	 * The minimum number of observations in a fit set is (the 2 is just for safety):
	 *	nrows - outerSegmentLength - innerSegmentLength - 2
	 */
	arma::uword minFitSetSize = this->nrows - outerSegmentLength - innerSegmentLength - 2;
	if(minFitSetSize <= this->maxNComp || this->maxNComp == 0) {
		this->maxNComp = minFitSetSize;
	}

	arma::uword j, orem, olen = 0, irem, ilenFixed,	ilen;

	IF_DEBUG(GAout << "Initialize segments with a test set size of " << testSetSize << std::endl);

	for(uint16_t rep = 0; rep < this->numReplications; ++rep) {
		arma::uvec shuffledRowNumbers = rowNumbers.shuffleAll(rng);

		orem = outerSegmentLengthRem;
		if(orem > 0) {
			if(innerSegmentLengthRem == 0) {
				irem = this->innerSegments - 1;
			} else {
				irem = innerSegmentLengthRem - 1;
			}
		} else {
			irem = innerSegmentLengthRem;
		}

		for(i = 0; i < this->outerSegments; ++i) {
			ilenFixed = innerSegmentLength;
			if(outerSegmentLength > 0) {
				if(orem > 0) {
					olen = outerSegmentLength + 1;
					ilenFixed = innerSegmentLengthBigOuter;
					--orem;
				} else {
					olen = outerSegmentLength;
					irem = innerSegmentLengthRem;
				}
			}

			n = 0;

			for(j = 0; j < this->innerSegments; ++j) {
				if(j < irem) {
					ilen = ilenFixed + 1;
				} else {
					ilen = ilenFixed;
				}

				/*
				 *
				 */
				arma::uvec inSegment(this->nrows - olen - ilen);

				if(n > 0) {
					inSegment.rows(0, n - 1) = shuffledRowNumbers.rows(0, n - 1);
				}

				if(n < this->nrows - olen - ilen) {
					inSegment.rows(n, inSegment.n_elem - 1) = shuffledRowNumbers.rows(n + ilen, this->nrows - olen - 1);
				}

				std::sort(inSegment.begin(), inSegment.end());

				/* First push back training set */
				this->segmentation.push_back(inSegment);
				IF_DEBUG(this->segmentation.back().t().raw_print(GAout, "Training set:"));
				/* Then add test set */
				this->segmentation.push_back(arma::sort(shuffledRowNumbers.rows(n, n + ilen - 1)));
				IF_DEBUG(this->segmentation.back().t().raw_print(GAout, "Test set:"));

				n += ilen;
			}

			/*
			 * "outer" segment
			 */
			if(olen > 0) {
				/* First push back training set */
				this->segmentation.push_back(arma::sort(shuffledRowNumbers.rows(0, n - 1)));
				IF_DEBUG(this->segmentation.back().t().raw_print(GAout, "Outer training set:"));
				/* Then add test set */
				this->segmentation.push_back(arma::sort(shuffledRowNumbers.rows(n, this->nrows - 1)));
				IF_DEBUG(this->segmentation.back().t().raw_print(GAout, "Outer test set:"));
			}

			/*
			 * Rotate shuffled row numbers (put the last segment in front)
			 */
			shuffledRowNumbers = arma::join_cols(shuffledRowNumbers.rows(n, this->nrows - 1), shuffledRowNumbers.rows(0, n - 1));
		}
	}
}

/**
 * Estimate the SEP for the given row ordering
 */
double PLSEvaluator::estSEP(uint16_t maxNComp) {
	double sumSEP = 0.0;
	OnlineStddev fitRSS(maxNComp);

	double cutoff;
	arma::uword optNComp;

	arma::vec residuals;
	arma::mat leftOutX;
	arma::vec leftOutY;
	OnlineStddev predSD;

	uint16_t rep = 0, outer, seg, comp;
	std::vector<arma::uvec>::const_iterator segmentIter = this->segmentation.begin();
	arma::uword j;

	while(rep++ < this->numReplications) {
		outer = 0;
		predSD.reset();

		while(outer++ < this->outerSegments) {
			/* Reset the variables to 0 */
			seg = 0;
			fitRSS.reset();

			/*
			 * Fit PLS models to predict the values in each segment once
			 */
			while(seg++ < this->innerSegments) {
				/* Segmentation iterator currently points to the `fit` segment */
				this->pls->viewSelectRows(*(segmentIter));
				this->pls->fit(maxNComp);

				/* Increment segmentation iterator to point to the `predict` segment */
				++segmentIter;

				leftOutY = this->pls->getY().rows(*segmentIter);
				leftOutX = this->pls->getXColumnView().rows(*segmentIter);

				for(comp = 0; comp < maxNComp; ++comp) {
					residuals = leftOutY - this->pls->predict(leftOutX, comp + 1);
					for(j = 0; j < residuals.n_elem; ++j) {
						fitRSS.update(residuals[j] * residuals[j], comp);
					}
				}

				/* Increment segmentation iterator to point to the next `fit` segment */
				++segmentIter;
			}

			/*
			 * Find best number of components based on the RSS plus one standard deviation
			 */
			IF_DEBUG(
				GAout << "EVALUATOR: MSE and SD" << std::endl;
				for(j = 0; j < maxNComp; ++j) {
					GAout << "\t (" << j + 1 << " comps.): " << fitRSS.mean(j) << " +- " << fitRSS.stddev(j) << std::endl;
				}
			)

			cutoff = fitRSS.mean(0);
			optNComp = 0;
			
			for(comp = 1; comp < maxNComp; ++comp) {
				if(fitRSS.mean(comp) < cutoff) {
					optNComp = comp;
					cutoff = fitRSS.mean(comp);
				}
			}

			IF_DEBUG(GAout << "EVALUATOR: Nr. of components with min. MSE: " << optNComp + 1 << " (max. " << maxNComp << ")" << std::endl)
			
			cutoff += fitRSS.stddev(optNComp);

			/* Go backwards again until the MSE is larger than the cutoff value */
			if(optNComp == 0) {
				optNComp = 1;
			} else {
				while(optNComp > 0 && fitRSS.mean(optNComp - 1) <= cutoff) {
					--optNComp;
				}
				++optNComp;
			}

			IF_DEBUG(GAout << "EVALUATOR: Opt. num. of components: " << optNComp << " (max. " << maxNComp << ")" << std::endl)

			/*
			 * Predict last segment with a model fit to the other observations using optNComp components
			 * (segmentation iterator points to the `test fit` segment
			 */
			this->pls->viewSelectRows(*segmentIter);
			this->pls->fit(optNComp);


			/* Increment the segmentation iterator to point to the `test predict` segment */
			++segmentIter;

			leftOutX = this->pls->getXColumnView().rows(*segmentIter);
			leftOutY = this->pls->getY().rows(*segmentIter);
			residuals = leftOutY - this->pls->predict(leftOutX, optNComp);

			/* Increment the segmentation iterator to point to the next `fit` segment */
			++segmentIter;

			predSD.update(residuals);
		}

		/*
		 * Calculate standard deviation of the residuals
		 */
		IF_DEBUG(GAout << "EVALUATOR: Resulting SEP: " << predSD.stddev() << std::endl)
		sumSEP -= predSD.stddev();
	}

	return sumSEP;
}

Evaluator* PLSEvaluator::clone() const {
	return new PLSEvaluator(*this);
}



