//
//  SynchronizedUnifGenerator.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 14.07.2013.
//
//
#include "config.h"

#ifdef HAVE_PTHREAD_H
#include "SynchronizedUnifGenerator__0__1.h"
#include <pthread.h>

pthread_mutex_t SynchronizedUnifGenerator__0__1::runifMutex = PTHREAD_MUTEX_INITIALIZER;
Rcpp::stats::UnifGenerator__0__1 SynchronizedUnifGenerator__0__1::unifGen;

double SynchronizedUnifGenerator__0__1::operator()() const {
	double rand = 0.0;
	int pthreadRC = pthread_mutex_lock(&SynchronizedUnifGenerator__0__1::runifMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)

	rand = SynchronizedUnifGenerator__0__1::unifGen();

	pthreadRC = pthread_mutex_unlock(&SynchronizedUnifGenerator__0__1::runifMutex);
	CHECK_PTHREAD_RETURN_CODE(pthreadRC)
	
	return rand;
}
#endif
