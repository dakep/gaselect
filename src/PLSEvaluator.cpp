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

	this->initSegmentation(testSetSize, seed);
}

PLSEvaluator::PLSEvaluator(const PLSEvaluator &other) :
	Evaluator(other.verbosity), numReplications(other.numReplications), outerSegments(other.outerSegments),
	innerSegments(other.innerSegments), nrows(other.nrows), cloned(true), segmentation(other.segmentation),
	maxNComp(other.maxNComp), minSegmentLength(other.minSegmentLength)
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
	ShuffledSet rowNumbers(nrows);
	arma::uword i, n = 0;

	if(testSetSize == 0.0 || this->outerSegments == 1) {
		testSetSize = 1.0 / (this->innerSegments + 1.0);
	}

	/*
	 * The size of the outer segment and the number of outer segments with one extra
	 * observation
	 */
	arma::uword outerSegmentLength = nrows * testSetSize;
	arma::uword outerSegmentLengthRem = nrows % outerSegmentLength;

	this->segmentation.reserve(2 * this->numReplications * (this->innerSegments + 1) * this->outerSegments);

	/*
	 * The size of the inner segment and the number of inner segments with one extra
	 * observation, if the smaller outer segment is used.
	 * If the larger outer segment is used, one inner segment with an extra observation
	 * won't get this extra observation. If no inner segment has an extra observation,
	 * one inner segment will have one observation less.
	 */
	arma::uword innerSegmentLength = (nrows - outerSegmentLength) / this->innerSegments;
	arma::uword innerSegmentLengthRem = (nrows - outerSegmentLength) % this->innerSegments;
	arma::uword innerSegmentLengthBigOuter = innerSegmentLength;

	if(outerSegmentLengthRem > 0 && innerSegmentLengthRem == 0) {
		--innerSegmentLengthBigOuter;
	}

	/*
	 * Update the maximal number of components, as innerSegmentLengthBigOuter is the minimal
	 * segment length.
	 */
	if(innerSegmentLengthBigOuter <= this->maxNComp) {
		this->minSegmentLength = innerSegmentLengthBigOuter;
		this->maxNComp = innerSegmentLengthBigOuter - 1;
	}

	arma::uword j, orem, olen = 0, irem, ilenFixed,	ilen;

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
				arma::uvec inSegment(nrows - olen - ilen);

				if(n > 0) {
					inSegment.rows(0, n - 1) = shuffledRowNumbers.rows(0, n - 1);
				}

				if(n < nrows - olen - ilen) {
					inSegment.rows(n, inSegment.n_elem - 1) = shuffledRowNumbers.rows(n + ilen, nrows - olen - 1);
				}

				std::sort(inSegment.begin(), inSegment.end());

				/* First push back training set */
				this->segmentation.push_back(inSegment);
				IF_DEBUG(this->segmentation.back().t().raw_print("Training set:"));
				/* Then add test set */
				this->segmentation.push_back(arma::sort(shuffledRowNumbers.rows(n, n + ilen - 1)));
				IF_DEBUG(this->segmentation.back().t().raw_print("Test set:"));

				n += ilen;
			}

			/*
			 * "outer" segment
			 */
			if(olen > 0) {
				/* First push back training set */
				this->segmentation.push_back(arma::sort(shuffledRowNumbers.rows(0, n - 1)));
				IF_DEBUG(this->segmentation.back().t().raw_print("Outer training set:"));
				/* Then add test set */
				this->segmentation.push_back(arma::sort(shuffledRowNumbers.rows(n, nrows - 1)));
				IF_DEBUG(this->segmentation.back().t().raw_print("Outer test set:"));
			}

			/*
			 * Rotate shuffled row numbers (put the last segment in front)
			 */
			shuffledRowNumbers = arma::join_cols(shuffledRowNumbers.rows(n, nrows - 1), shuffledRowNumbers.rows(0, n - 1));
		}
	}
}

/**
 * Estimate the SEP for the given row ordering
 */
double PLSEvaluator::estSEP(uint16_t maxNComp) {
	double sumSEP = 0.0;
	// (online) Sum of squares of differences from the (current) mean (residMeans)
	// M_2,n = sum( (x_i - mean_n) ^ 2 )
	arma::vec RSS(maxNComp);
	double RSSSDdelta = 0.0;
	arma::vec RSSSDm2n(maxNComp);
	std::vector<double> RSSSDmean(maxNComp, 0.0);
	double r2;

	double cutoff;
	arma::uword optNComp;

	arma::mat residuals;
	arma::mat leftOutX;
	arma::mat leftOutY;
	
	RSSSDm2n.zeros();
	RSS.zeros();

	uint16_t rep = 0, outer, seg, comp;
	arma::uword i = 0, j, n;

	while(rep++ < this->numReplications) {
		outer = 0;
		while(outer++ < this->outerSegments) {
			seg = 0;
			n = 0;

			/*
			 * Fit PLS models to predict the values in each segment once
			 */
			while(seg++ < this->innerSegments) {
				this->pls->viewSelectRows(this->segmentation[i++]);
				this->pls->fit(maxNComp);

				leftOutY = this->pls->getY().rows(this->segmentation[i]);
				leftOutX = this->pls->getXColumnView().rows(this->segmentation[i]);

				for(comp = 1; comp <= maxNComp; ++comp) {
					residuals = leftOutY - this->pls->predict(leftOutX, comp);

					for(j = 0; j < residuals.n_elem; ++j) {
						r2 = residuals[j] * residuals[j];
						RSS[comp] += r2;
						
						RSSSDdelta = r2 - RSSSDmean[comp];
						RSSSDmean[comp] = RSSSDmean[comp] + RSSSDdelta / (n + 1 + j);
						RSSSDm2n[comp] = RSSSDm2n[comp] + RSSSDdelta * (r2 - RSSSDmean[comp]);
					}
				}

				++i;
				n += j;
			}


			/*
			 * Find best number of components based on the RSS plus one standard deviation
			 */
			cutoff = RSS[0];
			optNComp = 0;
			
			for(comp = 1; comp < maxNComp; ++comp) {
				if(RSS[comp] < cutoff) {
					optNComp = comp;
					cutoff = RSS[comp];
				}
			}
			
			cutoff += sqrt(RSSSDm2n[optNComp] / (this->nrows - this->minSegmentLength));
			
			optNComp = 0;
			while(RSS[optNComp] > cutoff) {
				++optNComp;
				if(optNComp > maxNComp) {
					optNComp = 0;
					break;
				}
			}
			++optNComp;

			/*
			 * @TODO Predict last segment with a model fit to the other observations using optNComp components
			 */
			this->pls->viewSelectRows(this->segmentation[i++]);
			this->pls->fit(optNComp);

			leftOutX = this->pls->getXColumnView().rows(this->segmentation[i]);
			leftOutY = this->pls->getY().rows(this->segmentation[i]);
			residuals = leftOutY - this->pls->predict(leftOutX, optNComp);

			IF_DEBUG(GAout << "EVALUATOR: Resulting SEP:" << std::endl << arma::stddev(residuals.col(0)) << std::endl)

			sumSEP += arma::stddev(residuals.col(0));
		}
	}

	return -sumSEP;
}

Evaluator* PLSEvaluator::clone() const {
	return new PLSEvaluator(*this);
}



