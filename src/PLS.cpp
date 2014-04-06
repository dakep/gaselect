//
//  PLS.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 15.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "config.h"
#include "Logger.h"
#include "PLS.h"

void PLS::viewSelectColumns(const arma::uvec &columns) {
	this->viewXCol = this->X.cols(columns);

	this->currentViewState = COLUMNS;
}

void PLS::viewSelectRows(const arma::uvec &rows) {
	this->viewX = this->viewXCol.rows(rows);
	this->viewY = this->Y.rows(rows);
	this->currentViewState = ROWS;
}

arma::mat PLS::predict(const arma::mat &newX, uint16_t ncomp) const {
	const arma::cube& coefs = this->getCoefficients();
	const arma::mat& intercepts = this->getIntercepts();
	if(ncomp > coefs.n_slices) {
		GAerr << "Trying to predict with " << ncomp << " components when only " << coefs.n_slices << " components are available" << std::endl;
		throw Rcpp::exception("Can not predict values for a model with more components than fit components", __FILE__, __LINE__);
	}
	
	--ncomp;
	arma::mat pred = newX * coefs.slice(ncomp);
	pred.each_row() += intercepts.row(ncomp);
	return pred;
}

arma::cube PLS::predict(const arma::mat &newX) const {
	const arma::cube& coefs = this->getCoefficients();
	const arma::mat& intercepts = this->getIntercepts();
	arma::cube pred(newX.n_rows, coefs.n_cols, coefs.n_slices);
	for(uint16_t i = 0; i < coefs.n_slices; ++i) {
		pred.slice(i) = newX * coefs.slice(i);
		pred.slice(i).each_row() += intercepts.row(i);
	}

	return pred;
}

PLS* PLS::getInstance(PLSMethod method, const arma::mat &X, const arma::mat &Y) {
	PLS *ret;
	switch(method) {
		case SIMPLS:
			ret = new PLSSimpls(X, Y);
			break;
	}

	return ret;
}


