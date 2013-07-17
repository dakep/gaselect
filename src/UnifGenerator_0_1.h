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

class UnifGenerator_0_1 {
public:
	UnifGenerator_0_1() {};
	virtual ~UnifGenerator_0_1() {};
	
	inline virtual double operator()() {
		return UnifGenerator_0_1::unifGen();
	}
protected:
	static Rcpp::stats::UnifGenerator__0__1 unifGen;
};


#ifdef HAVE_PTHREAD_H
class SynchronizedUnifGenerator_0_1 : public UnifGenerator_0_1 {
public:
	SynchronizedUnifGenerator_0_1(const uint16_t bufferSize = 0) : bufferSize(bufferSize), curRandPos(bufferSize) {
		this->randBuffer.reserve(this->bufferSize);
	};
	~SynchronizedUnifGenerator_0_1() {}
	
	double operator()();
	
	SynchronizedUnifGenerator_0_1& operator=(const SynchronizedUnifGenerator_0_1& ) {
		throw std::logic_error("The synchronized uniform generator can not be assigned");
	}
	
private:
	static pthread_mutex_t runifMutex;

	const uint16_t bufferSize;
	uint16_t curRandPos;
	std::vector<double> randBuffer;
};
#endif

#endif // header def