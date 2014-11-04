//
//  PLSSimpls.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "config.h"

#include <iostream>
#include <stdexcept>
#include <RcppArmadillo.h>
#include "PLSSimpls.h"

const double PLSSimpls::NORM_TOL = 1e-20;

PLSSimpls::PLSSimpls(const arma::mat &X, const arma::vec &Y) : PLS(X, Y) {
	if(Y.n_cols > 1) {
		throw std::invalid_argument("The size of the seed must not be smaller than the RNG's seed size");
	}
}

PLSSimpls::~PLSSimpls() {
}

PLS* PLSSimpls::clone() const {
	PLSSimpls* clone = new PLSSimpls(arma::mat(this->X), arma::mat(this->Y));
	
//	clone->resultNComp = this->resultNComp;
//	clone->fittedValues = this->fittedValues;
//	clone->coef = this->coef;
//	clone->intercepts = this->intercepts;
//	clone->Ymean = this->Ymean;
//	clone->Xmean = this->Xmean;
	
	return clone;
}

inline void PLSSimpls::centerView() {
	// center X and Y

	switch(this->currentViewState) {
		case PLS::UNKNOWN:
			this->viewX = this->X;
			this->viewY = this->Y;
			break;
		case PLS::COLUMNS:
			this->viewX = this->viewXCol;
			this->viewY = this->Y;
			break;
		default:
			break;
	}

	this->Xmean = arma::mean(this->viewX);
	this->Ymean = arma::mean(this->viewY);

	this->viewX.each_row() -= this->Xmean;
	this->viewY -= this->Ymean;
}

/**
 * ncomp is 1-based!!!
 */
void PLSSimpls::fit(uint16_t ncomp) {
	uint16_t maxNComp = ((this->viewX.n_cols < this->viewX.n_rows) ? this->viewX.n_cols : this->viewX.n_rows - 1);
	if(ncomp == 0 || ncomp > maxNComp) {
		ncomp = maxNComp;
	}

	/*
	 * Init neccessary matrices and vectors
	 * Variable names are according to the original paper by S. de Jong (1993)
	 */

	this->coef.zeros(this->viewX.n_cols, ncomp);
	this->intercepts.zeros(ncomp);

	this->R.zeros(this->viewX.n_cols, ncomp);
	this->V.zeros(this->viewX.n_cols, ncomp);
	this->Qvec.zeros(ncomp);

	/*
	 * Center X and Y views
	 */
	this->centerView();

	arma::vec S = this->viewX.t() * this->viewY; // Cross product

	// Working vectors
	arma::vec t; // X block factor scores
	double tnorm = 1.0;
	arma::vec p; // X block factor loadings

	for(uint16_t i = 0; i < ncomp; ++i) {
		arma::vec r = this->R.unsafe_col(i); // X block factor weights
		arma::vec v = this->V.unsafe_col(i); // Orthogonal loadings

		t = this->viewX * S;

		t = t - arma::mean(t); // Center y block factor scores
		tnorm = arma::norm(t, 2); // Calculate norm

		// The norm of t can be zero (or close to it). This is unacceptable.
		if (tnorm < PLSSimpls::NORM_TOL) {
			throw std::underflow_error("All block-factor scores are (almost) zero.");
		}

		t = t / tnorm;  // Normalize scores
		r = S / tnorm;

		p = this->viewX.t() * t; // Calculate x loadings
		this->Qvec[i] = arma::dot(this->viewY, t); // Calculate y loadings

		if(i > 0) {
			v = arma::normalise(p - (this->V * this->V.t() * p)); // Make v orthogonal to previous loadings
		} else {
			v = arma::normalise(p);
		}

		S -= v * v.t() * S; // deflate S

		// R and V must not be updated because r resp. v are already pointers to the correct column

		this->coef.col(i) = this->R.cols(0, i) * this->Qvec.rows(0, i);
		this->intercepts[i] = this->Ymean - arma::dot(this->Xmean, this->coef.col(i));
	}

	this->resultNComp = ncomp;
}
