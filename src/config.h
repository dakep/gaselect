//
//  config.h
//  GenAlgTest
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


#ifdef HAVE_UINT64_MAX
	typedef uint64_t IntChromosome;
	#define INT_CHROMOSOME_MAX UINT64_MAX
#elif defined HAVE_UINT32_MAX
	typedef uint32_t IntChromosome;
	#define INT_CHROMOSOME_MAX UINT32_MAX
#else
#error "No suitable integral type supported"
#endif


// Mutation algorithm
// If the ratio of set to unset bits is greater than this number, a random position
// is likely to result in a bit with the desired state
#define RATIO_RANDOM_SEARCH 0.2

#endif
