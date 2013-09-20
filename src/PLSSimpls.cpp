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

PLSSimpls::PLSSimpls(const arma::mat &X, const arma::mat &Y, const bool fitValues) : PLS(X, Y, fitValues) {
	if(Y.n_cols > 1) {
		throw std::invalid_argument("The size of the seed must not be smaller than the RNG's seed size");
	}
	this->subviewChanged();
}

PLSSimpls::~PLSSimpls() {
}

PLS* PLSSimpls::clone() const {
	PLSSimpls* clone = new PLSSimpls(arma::mat(this->X), arma::mat(this->Y), this->fitValues);
	
	clone->viewX = this->viewX;
	clone->viewY = this->viewY;
	clone->viewXCol = this->viewXCol;
	
//	clone->resultNComp = this->resultNComp;
//	clone->fittedValues = this->fittedValues;
//	clone->coef = this->coef;
//	clone->intercepts = this->intercepts;
//	clone->Ymean = this->Ymean;
//	clone->Xmean = this->Xmean;
	
	return clone;
}

void PLSSimpls::subviewChanged() {
	PLS::subviewChanged();

	// center X and Y
	this->Xmean = ((arma::mat) arma::mean(this->viewX)).row(0);
	this->Ymean = ((arma::mat) arma::mean(this->viewY)).row(0);

	this->viewX.each_row() -= this->Xmean;
	this->viewY.each_row() -= this->Ymean;
}

void PLSSimpls::fit(uint16_t ncomp) {
	uint16_t maxNComp = ((this->viewX.n_cols < this->viewX.n_rows) ? this->viewX.n_cols : this->viewX.n_rows - 1);
	if(ncomp == 0 || ncomp > maxNComp) {
		ncomp = maxNComp;
	}

	// Init neccessary matrices and vectors
	// Variable names are accoridng to the original paper by S. de Jong (1993)

	this->coef.zeros(this->viewX.n_cols, this->viewY.n_cols, ncomp);
	this->intercepts.zeros(ncomp, this->viewY.n_cols);

	arma::mat TT; // X factor scores (only really needed if fitted values should be calculated)

	this->R.zeros(this->viewX.n_cols, ncomp);
	this->V.zeros(this->viewX.n_cols, ncomp);
	this->tQ.zeros(ncomp, this->viewY.n_cols);

	if(this->fitValues) {
		TT.zeros(this->viewX.n_rows, ncomp);
		this->fittedValues.zeros(this->viewY.n_rows, this->viewY.n_cols, ncomp);
	} else {
		this->fittedValues.zeros(1, 1, 1);
	}

	arma::mat S = this->viewX.t() * this->viewY; // Cross product

	// Working vectors
	arma::vec t; // X block factor scores
	double tnorm = 1.0;
	arma::vec p; // X block factor loadings
	arma::vec q; // Y block factor loadings / weights

	for(uint16_t i = 0; i < ncomp; ++i) {
		arma::vec r = this->R.unsafe_col(i); // X block factor weights
		arma::vec v = this->V.unsafe_col(i); // Orthogonal loadings

		r = S; // Only univariate responses are supported

		t = this->viewX * r;

		t = t - arma::mean(t); // Center y block factor scores
		tnorm = ( arma::sqrt(t.t() * t)[0] ); // Calculate norm
		t = t / tnorm;  // Normalize scores
		r = r / tnorm;

		p = this->viewX.t() * t; // Calculate x loadings
		q = this->viewY.t() * t; // Calculate y loadings

		if(i > 0) {
			v = p - (V * V.t() * p); // Make v orthogonal to previous loadings
		} else {
			v = p;
		}

		v = v / ( arma::sqrt(v.t() * v)[0] ); // Normalize orthogonal loadings

		S = S - v * v.t() * S; // deflate S
		this->tQ.row(i) = q.t();

		// R and V must not be updated because r resp. v are already pointers to the correct column

		this->coef.slice(i) = this->R.cols(0, i) * this->tQ.rows(0, i);
		this->intercepts.row(i) = this->Ymean - (this->Xmean * this->coef.slice(i));

		if(this->fitValues) {
			TT.col(i) = t;
			this->fittedValues.slice(i) = TT.cols(0, i) * tQ.rows(0, i);
			this->fittedValues.slice(i).each_row() += this->Ymean;
		}
	}

	this->resultNComp = ncomp;
	this->validResultState = true;
}