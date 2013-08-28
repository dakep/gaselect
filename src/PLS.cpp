//
//  PLS.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 15.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#include "config.h"
#include "PLS.h"

void PLS::subviewChanged() {
	this->validResultState = false;
}

void PLS::setSubmatrixView(const arma::uvec &rows, const arma::uvec &columns) {
	// As usually the submatrices are used alot and the columns and rows
	// are generally not consequtive, it is faster to create a copy of the
	// submatrix and work with the copy than to always access the submatrix view

	this->viewXCol = this->X.cols(columns);
	this->viewX = this->viewXCol.rows(rows);
	this->viewY = this->Y.rows(rows);

	this->subviewChanged();
}

void PLS::setSubmatrixViewColumns(const arma::uvec &columns) {
	// As usually the submatrices are used alot and the columns and rows
	// are generally not consequtive, it is faster to create a copy of the
	// submatrix and work with the copy than to always access the submatrix view

	this->viewX = this->viewXCol = this->X.cols(columns);
	this->viewY = this->Y;

	this->subviewChanged();
}

void PLS::setSubmatrixViewRows(const arma::uvec &rows, bool keepOldColumns) {
	// As usually the submatrices are used alot and the columns and rows
	// are generally not consequtive, it is faster to create a copy of the
	// submatrix and work with the copy than to always access the submatrix view

	if(keepOldColumns == true) {
		this->viewX = this->viewXCol.rows(rows);
	} else {
		this->viewX = this->X.rows(rows);
	}
	this->viewY = this->Y.rows(rows);

	this->subviewChanged();
}

// ncomp should be zero based
arma::mat PLS::predict(arma::mat newX, uint16_t ncomp) const {
	arma::cube coefs = this->getCoefficients();
	if(ncomp > coefs.n_slices) {
		Rcpp::Rcout << "Trying to predict with " << ncomp << " components when only " << coefs.n_slices << " components are available" << std::endl;
		throw Rcpp::exception("Can not predict values for a model with more components than fit components", __FILE__, __LINE__);
	}
	
	--ncomp;
	arma::mat pred = newX * coefs.slice(ncomp);
	pred.each_row() += this->getIntercepts().row(ncomp);
	return pred;
}

arma::cube PLS::predict(arma::mat newX) const {
	arma::cube coefs = this->getCoefficients();
	arma::mat intercepts = this->getIntercepts();
	arma::cube pred(newX.n_rows, coefs.n_cols, coefs.n_slices);
	for(uint16_t i = 0; i < coefs.n_slices; ++i) {
		pred.slice(i) = newX * coefs.slice(i);
		pred.slice(i).each_row() += intercepts.row(i);
	}

	return pred;
}

PLS* PLS::getInstance(PLSMethod method, const arma::mat &X, const arma::mat &Y, const bool fitValues) {
	PLS *ret;
	switch(method) {
		case SIMPLS:
			ret = new PLSSimpls(X, Y, fitValues);
			break;
	}

	return ret;
}


