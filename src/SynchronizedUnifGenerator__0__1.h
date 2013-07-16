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
#include <vector>
#endif

#include <RcppArmadillo.h>
#include <Rcpp/stats/random/runif.h>

#if RCPP_VERSION < Rcpp_Version(0, 10, 1)
namespace Rcpp {
	namespace stats {
		class UnifGenerator__0__1 : public ::Rcpp::Generator<false, double> {
		public:
			UnifGenerator__0__1( double _min = 0.0, double _max = 1.0) {}
			
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
	SynchronizedUnifGenerator__0__1(const uint16_t bufferSize = 0) : bufferSize(bufferSize), curRandPos(bufferSize) {
		this->randBuffer.reserve(this->bufferSize);
	};
	
	double operator()();
	
	SynchronizedUnifGenerator__0__1& operator=(const SynchronizedUnifGenerator__0__1& ) {
		throw std::logic_error("The synchronized uniform generator can not be assigned");
	}
	
private:
	static Rcpp::stats::UnifGenerator__0__1 unifGen;
	static pthread_mutex_t runifMutex;

	const uint16_t bufferSize;

	uint16_t curRandPos;
	std::vector<double> randBuffer;
};

#else
class SynchronizedUnifGenerator__0__1 {
public:
	SynchronizedUnifGenerator__0__1() {}

	inline double SynchronizedUnifGenerator__0__1::operator()() {
		return SynchronizedUnifGenerator__0__1::unifGen();
	}
private:
	static Rcpp::stats::UnifGenerator__0__1 unifGen;
};

#endif

#endif // header def