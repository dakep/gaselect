//
//  PLSSimpls.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 15.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef GenAlgPLS_PLSSimpls_h
#define GenAlgPLS_PLSSimpls_h

#include "config.h"

#include <RcppArmadillo.h>
#include "PLS.h"

class PLSSimpls : public PLS {
public:
	PLSSimpls(const arma::mat &X, const arma::mat &Y, const bool fitValues = true);
	~PLSSimpls();

	void fit(uint16_t ncomp = 0);
	arma::cube getCoefficients() const { return this->coef; }
	arma::mat getIntercepts() const { return this->intercepts; };
	arma::cube getFittedValues() const { return this->fittedValues; }

//	arma::mat getScores() { return this->TT; }
//	arma::mat getYLoadings() { return this->tQ.t(); }
//	arma::mat getProjection() { return this->R; }

	virtual PLS* clone() const;
	
protected:
	void subviewChanged();

private:
	arma::cube fittedValues;
	arma::cube coef; // Cube with coefficients
	arma::mat intercepts; // n x ncomp matrix with intercept terms
	arma::rowvec Ymean; // Column means of Y
	arma::rowvec Xmean; // Column means of X
};


#endif
