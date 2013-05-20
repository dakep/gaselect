//
//  PLSEvaluator.cpp
//  GenAlgTest
//
//  Created by David Kepplinger on 03.05.2013.
//
//

#include "config.h"

#include <algorithm>
#include "PLSEvaluator.h"

double PLSEvaluator::evaluate(Chromosome &ch) const {
	arma::uvec columnSubset = ch.toColumnSubset();
	// -2 because if the segmentLength would not be an exact integer some segments are one longer than others
	uint16_t maxNComp = ((columnSubset.n_elem < this->nrows) ? columnSubset.n_elem : this->nrows - this->segmentLength - 2);
	arma::vec sumSSD(maxNComp);
	arma::uvec rowNumbers = this->initRowNumbers();
	arma::uword rep = 0;
	sumSSD.zeros();
	
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->verbosity == DEBUG_VERBOSE) {
		Rcpp::Rcout << "EVALUATOR: Testing model with variables" << std::endl << columnSubset.t() << std::endl;
	}
#endif

	for(rep = 0; rep < this->numReplications; ++rep) {
		sumSSD += this->calcSSD(columnSubset, maxNComp, rowNumbers);
	}
	
	double bestFitness = -sumSSD.min();
	ch.setFitness(bestFitness);

#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->verbosity == DEBUG_VERBOSE) {
		Rcpp::Rcout << "EVALUATOR: Sum of sqrt(sum of squared differences) for every number of components:" << std::endl << sumSSD.t() << std::endl;
	}
#endif
	return bestFitness;
}

double PLSEvaluator::evaluate(arma::uvec &columnSubset) const {
	uint16_t maxNComp = ((columnSubset.n_elem < this->nrows) ? columnSubset.n_elem : this->nrows - this->segmentLength - 2);
	arma::vec sumSSD(maxNComp);
	arma::uvec rowNumbers = this->initRowNumbers();
	arma::uword rep = 0;
	sumSSD.zeros();
		
	for(rep = 0; rep < this->numReplications; ++rep) {
		sumSSD += this->calcSSD(columnSubset, maxNComp, rowNumbers);
	}
	
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->verbosity == DEBUG_VERBOSE) {
		Rcpp::Rcout << "EVALUATOR: Sum of sqrt(sum of squared differences) for every number of components:" << std::endl << sumSSD << std::endl;
	}
#endif

	return -sumSSD.min();
}

inline arma::uvec PLSEvaluator::initRowNumbers() const {
	arma::uword i = 0, j = 1;
	arma::uvec rowNumbers(this->nrows);
	
	for(; j < this->nrows; i += 2, j += 2) {
		rowNumbers[i] = i;
		rowNumbers[j] = j;
	}
	if(i < this->nrows) {
		rowNumbers[i] = i;
	}
	
	return rowNumbers;
}

arma::vec PLSEvaluator::calcSSD(arma::uvec &columnSubset, uint16_t ncomp, arma::uvec &rowNumbers) const {
	// (online) Sum of squares of differences from the (current) mean (residMeans)
	// M_2,n = sum( (x_i - mean_n) ^ 2 )
	arma::vec residM2n(ncomp);
	std::vector<double> residMeans(ncomp, 0.0);
	double delta = 0.0;
	uint16_t seg = 0, comp = 0;
	arma::uword randPos = 0, i = 0;
	arma::uword n = 0, nSeg = 0;

	// If not all segments are the same length, segmentLength is the longer one
	// (i.e. the incomplete segments will be one element shorter)
	arma::uword segmentLength = this->segmentLength;
	int32_t completeSegments = this->completeSegments;
	arma::uvec segment;
	arma::uvec notSegment;
		
	arma::mat residuals;
	arma::mat leftOutX;
	arma::mat leftOutY;
	
	residM2n.zeros();
	
	if(this->completeSegments > 0) {
		++segmentLength;
	}

	for(; seg < this->numSegments; ++seg) {
		// Determine segment
		// Randomly permute the last (numSegments - seg) * segmentLength elements
		// in rowNumbers and swap them with the first segmentLength elements
		// --> first segmentLength elements are the left out elements
		// all other elements are included elements
		if(seg == this->numSegments - 1) {
			segment = rowNumbers.rows(this->nrows - segmentLength, this->nrows - 1);
			notSegment = rowNumbers.rows(0, this->nrows - segmentLength - 1);
		} else {
			for(i = 0; i < segmentLength; ++i) {
				// find a random position in the back of the array
				// (first elements are already used elements or current ones)
				// The probability that the uniform generator returns *exactly* 1 is zero, but it might happen
				// anyway. Substracting a very small number may result in negative results, but the
				// integer is unsigned so it can not get smaller than 0
				randPos = n + i + (arma::uword) (this->unifGen() * (this->nrows - n - i));
								
				std::swap(rowNumbers[n + i], rowNumbers[randPos]);
				std::swap(rowNumbers[i], rowNumbers[n + i]);
			}
			segment = rowNumbers.rows(0, segmentLength - 1);
			notSegment = rowNumbers.rows(segmentLength, this->nrows - 1);
		}
		
#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->verbosity == DEBUG_VERBOSE) {
			Rcpp::Rcout << "EVALUATOR: " << seg << ". (not)segment:" << std::endl << "\t" << segment.t() << std::endl << "\t" << notSegment.t() << std::endl << std::endl;
		}
#endif
		leftOutX = this->pls->getX().submat(segment, columnSubset);
		leftOutY = this->pls->getY().rows(segment);
		this->pls->setSubmatrixView(notSegment, columnSubset);
		
		if(--completeSegments == 0) {
			--segmentLength;
		}
		
		
		this->pls->fit(ncomp);

		// Calculate the standard error for observations not present in this segment
		for(comp = 0; comp < ncomp; ++comp) {
			residuals = leftOutY - this->pls->predict(leftOutX, comp);

			// Only consider multiple PLS (not multivariate), i.e. only use first response vector
			for (nSeg = 0; nSeg < residuals.n_rows; ++nSeg) {
				delta = residuals[nSeg] - residMeans[comp];
				residMeans[comp] = residMeans[comp] + delta / (n + 1 + nSeg);
				residM2n[comp] = residM2n[comp] + delta * (residuals[nSeg] - residMeans[comp]);
			}
		}
		
		n += nSeg;
	}
	
	
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->verbosity == DEBUG_VERBOSE) {
		Rcpp::Rcout << "EVALUATOR: Resulting M2n:" << std::endl << arma::sqrt(residM2n) << std::endl;
	}
#endif
	
//	Altough the square root doesn't change the order of the values,
//	it may be necessary to accurately reflect the "distance" between two values
	return arma::sqrt(residM2n);
}


