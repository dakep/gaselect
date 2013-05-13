//
//  config.h
//  GenAlgTest
//
//  Created by David Kepplinger on 08.05.2013.
//
//

#ifndef GenAlgPLS_config_h
#define GenAlgPLS_config_h

#include <inttypes.h>

// #define ENABLE_DEBUG_VERBOSITY

//#define TIMING_BENCHMARK

#define MIN_MUTATION_PROBABILITY 0.000001 // 1e-6
#define TAB_DELIMITER "    "
#define PRECISION 8
#define WIDTH PRECISION + 5
#define DELIMITER_POSITION 4
#define BITS_PER_BYTE 8
#define RNG_MAX_BITS 32

typedef uint64_t IntChromosome;
#define INT_CHROMOSOME_BITS 64
#define INT_CHROMOSOME_MAX UINT64_MAX
#define RANDS_PER_INT_CHROMOSOME 2


#endif
