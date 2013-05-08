//
//  PLS.h
//

#ifndef GenAlgPLS_PLS_h
#define GenAlgPLS_PLS_h

#include "config.h"

#include <RcppArmadillo.h>
#include <inttypes.h>

enum PLSMethod {
	SIMPLS = 0
};

class PLS {

public:
	PLS(const arma::mat &X, const arma::mat &Y, const bool fitValues = true);
	virtual ~PLS();
	
	static PLS* getInstance(PLSMethod method, const arma::mat &X, const arma::mat &Y, const bool fitValues = true);
	
	virtual void setSubmatrixView(const arma::uvec &rows, const arma::uvec &columns);
	virtual void setSubmatrixViewColumns(const arma::uvec &columns);
	virtual void setSubmatrixViewRows(const arma::uvec &rows);

	/**
	 * Fit a PLS model to the data with the previously set view
	 * with up to ncomp components
	 */
	virtual void fit(uint16_t ncomp = 0) =0;
	
	/**
	 * Get the coefficients of the last fit (i.e. coefficients
	 * that are obtained with ncomp specified in the last call
	 * to PLS::fit.
	 */
	virtual arma::cube getCoefficients() const = 0;
	
	/**
	 * Returns the intercept term for every number of components
	 * i.e. ncomp x nresp matrix
	 */
	virtual arma::mat getIntercepts() const = 0;
	
	/**
	 * Check whether the results returned by the other methods (e.g. getCoefficients, getIntercepts, ...)
	 * are valid.
	 *
	 * Those results may be invalidated by changing the subview
	 */
	bool isResultValid() const { return this->validResultState; }

	arma::uword getNumberOfPredictorVariables() const { return this->viewX.n_cols; }
	arma::uword getNumberOfResponseVariables() const { return this->viewY.n_cols; }
	arma::uword getNumberOfObservations() const { return this->viewX.n_rows; }

	/**
	 * Returns the number components the last fit was performed with
	 */
	uint16_t getResultNComp() const { return this->resultNComp; }
	
	// ncomp should be zero based
	arma::mat predict(arma::mat newX, uint16_t ncomp) const;
	arma::cube predict(arma::mat newX) const;

	const arma::mat getX() const { return this->X; }
	const arma::mat getY() const { return this->Y; }
protected:
	const arma::mat X;
	const arma::mat Y;
	const bool fitValues;
	
	bool validResultState;
	uint16_t resultNComp;

	arma::mat viewX;
	arma::mat viewY;

	virtual void subviewChanged();
};

#include "PLSSimpls.h"

#endif
