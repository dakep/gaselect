//
//  Logger.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 14.10.2013.
//
//

#ifndef GenAlgPLS_Logger_h
#define GenAlgPLS_Logger_h

#include <iostream>
#include <streambuf>
#include <string>

template <bool ERROR_STREAM>
class LoggerStreamBuffer : public std::streambuf {
public:
	LoggerStreamBuffer();
	virtual ~LoggerStreamBuffer();

	void flushThreadSafeBuffer();

	void enableThreadSafety(bool threadSafe) {
		this->flushThreadSafeBuffer();
		this->threadSafe = threadSafe;
	}
protected:
	virtual std::streamsize xsputn(const char *s, std::streamsize n);
	virtual int overflow(int c = EOF);
	virtual int sync();

private:
	bool threadSafe;
	std::string tsBuffer;
#ifdef HAVE_PTHREAD_H
	pthread_mutex_t printMutex;
#endif
};


template <bool ERROR_STREAM>
class Logger : public std::ostream {
private:
	typedef LoggerStreamBuffer<ERROR_STREAM> Buffer;
	Buffer* buf;

public:
	Logger() : std::ostream(new Buffer()), buf(static_cast<Buffer*>(rdbuf())) {};
	~Logger() {
		if(this->buf != NULL) {
			delete this->buf;
			this->buf = NULL;
		}
	}

	void flushThreadSafeBuffer() {
		if(this->buf != NULL) {
			this->buf->flushThreadSafeBuffer();
		}

	}

	void enableThreadSafety(bool threadSafe = true) {
		if(this->buf != NULL) {
			this->buf->enableThreadSafety(threadSafe);
		}
	}
};

extern Logger<false> GAout;
extern Logger<true> GAerr;

#endif
