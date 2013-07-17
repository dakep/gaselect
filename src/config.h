//
//  config.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 08.05.2013.
//
//

#ifndef GenAlgPLS_config_h
#define GenAlgPLS_config_h

#include "autoconfig.h"

#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif


// #define ENABLE_DEBUG_VERBOSITY

//#define TIMING_BENCHMARK

#define MIN_MUTATION_PROBABILITY 0.000001 // 1e-6
#define TAB_DELIMITER "    "
#define PRECISION 8
#define WIDTH PRECISION + 5
#define DELIMITER_POSITION 4
#define BITS_PER_BYTE 8
#define RNG_MAX_BITS 32
#define INT_RNG_MAX 4294967296 // 2^RNG_MAX_BITS

#ifdef HAVE_UNSIGNED_LONG_LONG
	typedef unsigned long long IntChromosome;
#elif
	typedef unsigned long IntChromosome;
#endif

#if (defined HAVE_CLIMITS || defined HAVE_LIMITS_H) 
	#ifdef HAVE_UNSIGNED_LONG_LONG
		#define INT_CHROMOSOME_MAX_VAL ULLONG_MAX
	#else
		#define INT_CHROMOSOME_MAX_VAL ULONG_MAX
	#endif
#endif

// Mutation algorithm
// If the ratio of set to unset bits is greater than this number, a random position
// is likely to result in a bit with the desired state
#define RATIO_RANDOM_SEARCH 0.2

#ifndef _REENTRANT
#undef HAVE_PTHREAD_H
#endif


//#undef HAVE_PTHREAD_H
//#define HAVE_PTHREAD_H 1


/**
 * Define the random numbers buffer size of the different uniform number generator instances
 */
#define UNIF_GENERATOR_BUFFER_SIZE_MAIN 8192 // this takes 64kiB on the heap!
#define UNIF_GENERATOR_BUFFER_SIZE_THREAD 8192 // this takes 64kiB on the heap!

// Buffer 256 (i.e. 2kiB) of random numbers
#define RANDOM_NUMBER_BUFFER 256

#ifdef ENABLE_DEBUG_VERBOSITY
#define CHECK_PTHREAD_RETURN_CODE(rc) if((rc) != 0) { Rcpp::Rcerr << "Warning: Call to pthread function failed with error code " << (rc) << " in " << __FILE__ << ":" << __LINE__ << std::endl; }
#else
#define CHECK_PTHREAD_RETURN_CODE(rc)
#endif

#endif
