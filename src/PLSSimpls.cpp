//
//  PLSSimpls.cpp
//  gaselect
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
	 * Center X and Y views
	 */
	this->centerView();

	/*
	 * Init neccessary matrices and vectors
	 * Variable names are according to the original paper by S. de Jong (1993)
	 */

	this->coef.zeros(this->viewX.n_cols, ncomp);
	this->intercepts.zeros(ncomp);
	this->V.set_size(this->viewX.n_cols, ncomp);

	arma::vec S = this->viewX.t() * this->viewY; // Cross product

	// Working vectors
	arma::vec t; // X block factor scores
	double tnorm = 1.0;
	double q;

	for(uint16_t i = 0; i < ncomp; ++i) {
		arma::vec v = this->V.unsafe_col(i);
		t = this->viewX * S;

		t -= arma::mean(t); // Center y block factor scores
		tnorm = arma::norm(t, 2); // Calculate norm

		// The norm of t can be zero (or close to it). This is unacceptable.
		if (tnorm < PLSSimpls::NORM_TOL) {
			throw std::underflow_error("All block-factor scores are (almost) zero.");
		}

		t /= tnorm;  // Normalize scores

		v = this->viewX.t() * t; // Calculate x loadings
		q = arma::dot(this->viewY, t); // Calculate y loadings

		if(i > 0) {
			v = v - this->V.cols(0, i - 1) * this->V.cols(0, i - 1).t() * v; // Make v orthogonal to previous loadings
			this->coef.col(i) = this->coef.col(i - 1) + (S * q / tnorm);
		} else {
			this->coef.col(i) = S * q / tnorm;
		}

		v /= arma::norm(v, 2);

		S = S - v * v.t() * S; // deflate S

		this->intercepts[i] = this->Ymean - arma::dot(this->Xmean, this->coef.col(i));
	}

	this->resultNComp = ncomp;
}
