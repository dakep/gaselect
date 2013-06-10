//
//  UnifGenerator__0__1.h
//
//  Created by David Kepplinger on 10.06.2013.
//
//

#ifndef GenAlgTest_UnifGenerator__0__1_h
#define GenAlgTest_UnifGenerator__0__1_h

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

#endif
