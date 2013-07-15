//
//  SynchronizedUnifGenerator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 14.07.2013.
//
//

#ifndef GenAlgPLS_SynchronizedUnifGenerator_h
#define GenAlgPLS_SynchronizedUnifGenerator_h

#include "config.h"

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#endif

#include <RcppArmadillo.h>
#include <Rcpp/stats/random/runif.h>

#if RCPP_VERSION < Rcpp_Version(0, 10, 1)
namespace Rcpp {
	namespace stats {
		class UnifGenerator__0__1 : public ::Rcpp::Generator<false, double> {
		public:
			UnifGenerator__0__1( double min_ = 0.0, double max_ = 1.0) {}
			
			inline double operator()() const {
				double u;
				do {u = unif_rand();} while (u <= 0 || u >= 1);
				return u;
			}
		} ;
	}
}
#endif

#ifdef HAVE_PTHREAD_H
class SynchronizedUnifGenerator__0__1 {
public:
	SynchronizedUnifGenerator__0__1() {};
	double operator()() const;
private:
	static Rcpp::stats::UnifGenerator__0__1 unifGen;
	static pthread_mutex_t runifMutex;
};

#else
typedef Rcpp::stats::UnifGenerator__0__1 SynchronizedUnifGenerator__0__1;
#endif

#endif
