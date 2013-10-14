//
//  Logger.cpp
//  GenAlgPLS
//
//  Created by David Kepplinger on 14.10.2013.
//
//

#include "config.h"

#include "Logger.h"
#include <RcppArmadillo.h>
#include <exception>

#ifdef HAVE_PTHREAD_H

#ifdef ENABLE_DEBUG_VERBOSITY
#define CHECK_PTHREAD_RETURN_CODE(expr) {int rc = expr; if((rc) != 0) { Rcpp::Rcout << "Warning: Call to pthread function failed with error code " << (rc) << " in " << __FILE__ << ":" << __LINE__ << std::endl; }}

#else
#define CHECK_PTHREAD_RETURN_CODE(expr) {expr;}
#endif

template <>
LoggerStreamBuffer<false>::LoggerStreamBuffer() : threadSafe(false) {
	int pthreadRC = pthread_mutex_init(&this->printMutex, NULL);
	if(pthreadRC != 0) {
		throw std::runtime_error("Mutex to synchronize printing could not be initialized");
	}
}

template <>
LoggerStreamBuffer<true>::LoggerStreamBuffer() : threadSafe(false) {
	int pthreadRC = pthread_mutex_init(&this->printMutex, NULL);
	if(pthreadRC != 0) {
		throw std::runtime_error("Mutex to synchronize printing could not be initialized");
	}
}

template <>
LoggerStreamBuffer<false>::~LoggerStreamBuffer() {
	pthread_mutex_destroy(&this->printMutex);
}

template <>
LoggerStreamBuffer<true>::~LoggerStreamBuffer() {
	pthread_mutex_destroy(&this->printMutex);
}

template <>
inline std::streamsize LoggerStreamBuffer<false>::xsputn(const char *s, std::streamsize n) {
	if(this->threadSafe) {
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		this->tsBuffer.append(s, n);
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
	} else {
		Rprintf("%.*s", n, s);
	}
	return n;
}

template <>
inline std::streamsize LoggerStreamBuffer<true>::xsputn(const char *s, std::streamsize n) {
	if(this->threadSafe) {
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		this->tsBuffer.append(s, n);
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
	} else {
		REprintf("%.*s", n, s);
	}
	return n;
}

template <>
int LoggerStreamBuffer<false>::overflow(int c) {
	if(c != EOF) {
		if(this->threadSafe) {
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
			this->tsBuffer.append(1, (char) c);
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
		} else {
			Rprintf("%.1s", &c);
		}
	}
	return c;
}

template <>
int LoggerStreamBuffer<true>::overflow(int c) {
	if(c != EOF) {
		if(this->threadSafe) {
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
			this->tsBuffer.append(1, (char) c);
			CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
		} else {
			Rprintf("%.1s", &c);
		}
	}
	return c;
}

template <>
int LoggerStreamBuffer<false>::sync() {
	if(!this->threadSafe) {
		R_FlushConsole();
	}
	return 0;
}

template <>
int LoggerStreamBuffer<true>::sync() {
	if(!this->threadSafe) {
		R_FlushConsole();
	}
	return 0;
}

template <>
void LoggerStreamBuffer<false>::flushThreadSafeBuffer() {
	if(this->threadSafe) {
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		Rprintf("%.*s", this->tsBuffer.length(), this->tsBuffer.c_str());
		R_FlushConsole();
		this->tsBuffer.clear();
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
	}
}

template <>
void LoggerStreamBuffer<true>::flushThreadSafeBuffer() {
	if(this->threadSafe) {
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_lock(&this->printMutex))
		Rprintf("%.*s", this->tsBuffer.length(), this->tsBuffer.c_str());
		R_FlushConsole();
		this->tsBuffer.clear();
		CHECK_PTHREAD_RETURN_CODE(pthread_mutex_unlock(&this->printMutex))
	}
}

#else

template <> LoggerStreamBuffer<false>::LoggerStreamBuffer() : threadSafe(false) {}
template <> LoggerStreamBuffer<true>::LoggerStreamBuffer() : threadSafe(false) {}

template <> LoggerStreamBuffer<false>::~LoggerStreamBuffer() {}
template <> LoggerStreamBuffer<true>::~LoggerStreamBuffer() {}


template <>
inline std::streamsize LoggerStreamBuffer<false>::xsputn(const char *s, std::streamsize n) {
	Rprintf("%.*s", n, s);
	return n;
}

template <>
inline std::streamsize LoggerStreamBuffer<true>::xsputn(const char *s, std::streamsize n) {
	REprintf("%.*s", n, s);
	return n;
}

template <>
inline int LoggerStreamBuffer<false>::overflow(int c) {
	if(c != EOF) {
		Rprintf("%.1s", &c);
	}
	return c;
}

template <>
inline int LoggerStreamBuffer<true>::overflow(int c) {
	if(c != EOF) {
		REprintf("%.1s", &c);
	}
	return c;
}

template <>
int LoggerStreamBuffer<false>::sync() {
	R_FlushConsole();
	return 0;
}

template <>
int LoggerStreamBuffer<true>::sync() {
	R_FlushConsole();
	return 0;
}

template <>
void LoggerStreamBuffer<false>::flushThreadSafeBuffer() {
}

template <>
void LoggerStreamBuffer<true>::flushThreadSafeBuffer() {
}

#endif
Logger<false> GAout;
Logger<true> GAerr;