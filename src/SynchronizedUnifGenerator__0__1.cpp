//
//  SynchronizedUnifGenerator.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 14.07.2013.
//
//
#include "config.h"
#include "SynchronizedUnifGenerator__0__1.h"

Rcpp::stats::UnifGenerator__0__1 SynchronizedUnifGenerator__0__1::unifGen;

#ifdef HAVE_PTHREAD_H
#include <pthread.h>
#include <exception>
#include <vector>

pthread_mutex_t SynchronizedUnifGenerator__0__1::runifMutex = PTHREAD_MUTEX_INITIALIZER;

double SynchronizedUnifGenerator__0__1::operator()() {
	if(this->bufferSize > 0) {
		if(this->curRandPos >= this->bufferSize) {
			/*
			 * All buffered random numbers have been used -- create new ones
			 */
			int pthreadRC = pthread_mutex_lock(&SynchronizedUnifGenerator__0__1::runifMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
			
			for(uint16_t i = 0; i < this->bufferSize; ++i) {
				this->randBuffer[i] = SynchronizedUnifGenerator__0__1::unifGen();
			}
			
			pthreadRC = pthread_mutex_unlock(&SynchronizedUnifGenerator__0__1::runifMutex);
			CHECK_PTHREAD_RETURN_CODE(pthreadRC)
			
			this->curRandPos = 0;
		}

		return this->randBuffer[this->curRandPos++];
	} else {
		double rand = 0.0;
		int pthreadRC = pthread_mutex_lock(&SynchronizedUnifGenerator__0__1::runifMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		rand = SynchronizedUnifGenerator__0__1::unifGen();
		
		pthreadRC = pthread_mutex_unlock(&SynchronizedUnifGenerator__0__1::runifMutex);
		CHECK_PTHREAD_RETURN_CODE(pthreadRC)
		
		return rand;
	}
}
#endif
