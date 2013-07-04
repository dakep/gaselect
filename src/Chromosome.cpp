//
//  Chromosome.cpp
//
//

#include "config.h"

#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <RcppArmadillo.h>
#include <Rcpp/stats/random/rbinom.h>

#include "Chromosome.h"
#include "GenAlg.h"

#if defined HAVE_STRINGS_H
#include <strings.h>
//#elif defined FFS_VS2005
//#include <intrin.h>
//#pragma intrinsic(_BitScanForward)
#endif

#include <sys/time.h>

using namespace Rcpp;

Chromosome::Chromosome(const Control &ctrl, VariablePositionPopulation &varPosPop) : ctrl(ctrl), tgeom(1. - ctrl.mutationProbability), varPosPop(varPosPop) {
	// Determine the number of IntChromosome bit values that are
	// needed to represent all genes
	this->numParts = (uint16_t) this->ctrl.chromosomeSize / Chromosome::BITS_PER_PART;

	int rest = this->ctrl.chromosomeSize % Chromosome::BITS_PER_PART;

	if(rest > 0) {
		this->numParts += 1;
		this->unusedBits = Chromosome::BITS_PER_PART - rest;
	} else {
		this->unusedBits = 0;
	}

	this->fitness = 0.0;

	// Initialize chromosome parts randomly
	this->chromosomeParts.resize(this->numParts, 0);
	this->initChromosomeParts();
}

Chromosome::Chromosome(const Chromosome &other, bool copyChromosomeParts) : ctrl(other.ctrl), tgeom(other.tgeom), varPosPop(other.varPosPop) {
	this->copyFrom(other, copyChromosomeParts);
}

void Chromosome::copyFrom(const Chromosome &other, bool copyChromosomeParts) {
	this->fitness = other.fitness;
	this->numParts = other.numParts;
	this->unusedBits = other.unusedBits;
	
	// Copy chromosome parts
	if(copyChromosomeParts) {
		this->chromosomeParts = other.chromosomeParts;
	} else {
		this->chromosomeParts.resize(this->numParts, 0);
	}
}

//Chromosome::~Chromosome() {
//}

inline void Chromosome::initChromosomeParts() {
#ifdef TIMING_BENCHMARK
	timeval start, end;

	gettimeofday(&start, NULL);
#endif
	uint16_t bitsToSet = this->ctrl.minVariables + this->unifGen() * (this->ctrl.maxVariables - this->ctrl.minVariables);

#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << "Init chromosome with " << bitsToSet << " bits set" << std::endl;
		Rcout << "First part before initializing: " << this->chromosomeParts[0] << std::endl;
	}
#endif

	VariablePositionPopulation::const_iterator setPosIter = this->varPosPop.shuffle(bitsToSet, this->unusedBits);

	uint16_t part = 0;
	uint16_t offset = 0;
	
	for(; setPosIter != this->varPosPop.end(); ++setPosIter) {
		part = (*setPosIter) / Chromosome::BITS_PER_PART;
		offset = (*setPosIter) % Chromosome::BITS_PER_PART;
		this->chromosomeParts[part] |= (((IntChromosome) 1) << offset);
	}
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << "Initialized chromosome with " << this->popcount() << " bits set" << std::endl;
	}
#endif
#ifdef TIMING_BENCHMARK

	gettimeofday(&end, NULL);
	Rcout << "Init chromosome took " << (end.tv_sec * 1000.0 + (end.tv_usec / 1000.0)) - (start.tv_sec * 1000.0 + (start.tv_usec / 1000.0)) << " milliseconds" << std::endl;

#endif
}

std::vector<Chromosome> Chromosome::mateWith(const Chromosome &other) {
#ifdef TIMING_BENCHMARK
	timeval start, end;

	gettimeofday(&start, NULL);
#endif
	if(other.ctrl.chromosomeSize != this->ctrl.chromosomeSize) {
		throw InvalidCopulationException(__FILE__, __LINE__);
	}

	std::vector<Chromosome> children;

	Chromosome child1(*this, false);
	Chromosome child2(*this, false);

	IntChromosome randomMask = 0;
	IntChromosome negRandomMask = 0;

	for(uint16_t i = 0; i < this->numParts; ++i) {
		if(this->chromosomeParts[i] == other.chromosomeParts[i]) {
			// Just copy the chromosome part to both children if it is the same
			// for both parents
			child1.chromosomeParts[i] = child2.chromosomeParts[i] = this->chromosomeParts[i];
			if(this->ctrl.verbosity >= MORE_VERBOSE) {
				Rcout << "Chromosome part is the same for both parents -- copy part to both children" << std::endl;
			}
		} else {
			// Randomly pick some bits from one chromosome and some bits from the other
			// chromosome
			do {
				randomMask = this->runif();
				if(i == 0) {
					randomMask <<= this->unusedBits;
				}
			} while(randomMask == 0);
			negRandomMask = ~randomMask;

			child1.chromosomeParts[i] = (this->chromosomeParts[i] & randomMask) | (other.chromosomeParts[i] & negRandomMask);
			child2.chromosomeParts[i] = (this->chromosomeParts[i] & negRandomMask) | (other.chromosomeParts[i] & randomMask);

#ifdef ENABLE_DEBUG_VERBOSITY
			if(this->ctrl.verbosity == DEBUG_VERBOSE) {
				Rcout << "Mask for part " << i << ": ";
				this->printBits(Rcout, randomMask, (i == 0) ? this->unusedBits : 0) << std::endl;

				Rcout << "Resulting parts:" << std::endl
				<< "Child 1:";
				this->printBits(Rcout, child1.chromosomeParts[i], (i == 0) ? this->unusedBits : 0) << std::endl
				<< "Child 2:";
				this->printBits(Rcout, child2.chromosomeParts[i], (i == 0) ? this->unusedBits : 0) << std::endl;
			}
#endif
		}
	}

	children.reserve(2);
	children.push_back(child1);
	children.push_back(child2);

#ifdef TIMING_BENCHMARK
	gettimeofday(&end, NULL);
	Rcout << "Copulating took " << (end.tv_sec * 1000.0 + (end.tv_usec / 1000.0)) - (start.tv_sec * 1000.0 + (start.tv_usec / 1000.0)) << " milliseconds" << std::endl;
#endif

	return children;
}

bool Chromosome::mutate() {
	static std::vector<uint16_t> positionPopulation(this->ctrl.chromosomeSize);

	uint16_t currentlySetBits = this->popcount();
	uint16_t currentlyUnsetBits = this->ctrl.chromosomeSize - currentlySetBits;

	int32_t numChangeBits = 0; // a negative value means unsetting bits and a positive value means setting new bits

	if((this->ctrl.minVariables - currentlySetBits) > 0) { // Too few bits are set -- we MUST add some
		numChangeBits = (this->ctrl.minVariables - currentlySetBits);		
	} else if((this->ctrl.maxVariables - currentlySetBits) < 0) { // Too many bits are set -- we MUST unset some
		numChangeBits = (this->ctrl.maxVariables - currentlySetBits);
	}

	// Now set a random number of bits and unset a random number of bits
	numChangeBits += this->tgeom(this->ctrl.maxVariables - currentlySetBits) - this->tgeom(currentlySetBits - this->ctrl.minVariables);

	if(this->ctrl.maxVariables - currentlySetBits > 0) {
		/*
		 * We may set some additional bits
		 * numChangeBits is 0 if the number of currently set bits is in the predefined range
		 * or *positive* if too few bits are set
		 */
		numChangeBits += this->tgeom(this->ctrl.maxVariables - currentlySetBits - numChangeBits);
	}

	if(currentlySetBits - this->ctrl.minVariables > 0) {
		/*
		 * We may unset some bits
		 * numChangeBits is 0 if the number of currently set bits is in the predefined range
		 * or *negative* if too many bits are set
		 */
		numChangeBits -= this->tgeom(currentlySetBits - this->ctrl.minVariables + numChangeBits);
	}

#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		if(numChangeBits != 0) {
			Rcout << "###########################" << std::endl
			<< "Changing " << numChangeBits << " bits" << std::endl
			<< "###########################" << std::endl
			<< "Before mutation:" << *this << std::endl;
		} else {
			Rcout << "Changing " << numChangeBits << " bits" << std::endl;
		}
		
	}
#endif

	if(numChangeBits == 0) {
		return false;
	} else if(numChangeBits < 0) {
		/*
		 * We have to unset -numChangeBits bits, i.e. flip from 1 to 0
		 */
		this->shuffle(positionPopulation, currentlySetBits, -numChangeBits);
		std::vector<uint16_t> removeBitsPos(positionPopulation.begin(), positionPopulation.begin() - numChangeBits);
		std::sort(removeBitsPos.begin(), removeBitsPos.end());
		std::vector<uint16_t>::iterator removeBitsPosIt = removeBitsPos.begin();

		IntChromosome mask = 0;
		IntChromosome tmp;
		int8_t totalShift = this->unusedBits - 1; // 0 based -- position 0 is the least significant bit
		int8_t trailingZeros = totalShift; // 1 based -- the value 1 means 1 trailing zero
		
		uint16_t onesCount = 0;
		
		for(uint16_t i = 0; (i < this->numParts) && (removeBitsPosIt != removeBitsPos.end()); ++i) {
			tmp = this->chromosomeParts[i];
			do {
				tmp >>= trailingZeros + 1;
				
				trailingZeros = this->ctz(tmp); // get number of 0's before the first 1 (starting with least significant bit)
				
				if(trailingZeros >= Chromosome::BITS_PER_PART - totalShift - 1) {
					trailingZeros = Chromosome::BITS_PER_PART - totalShift - 1;
				}
				
				totalShift += trailingZeros + 1;
				
				if(totalShift < Chromosome::BITS_PER_PART) {
					if(onesCount == (*removeBitsPosIt)) {
						//unset bit
						mask |= ((IntChromosome) 1) << totalShift;
						++removeBitsPosIt;
					}
					++onesCount; // increment after comparision because the position is 0 based
				}
			} while(totalShift < Chromosome::BITS_PER_PART);
			
			this->chromosomeParts[i] ^= mask; // toggle all bits set in the mask
			
			totalShift = -1;
			trailingZeros = -1;
			mask = 0;
		}
	} else { // (numChangeBits > 0)
		/*
		 * We have to set numChangeBits bits, i.e. flip from 0 to 1
		 */
		this->shuffle(positionPopulation, currentlyUnsetBits, numChangeBits);
		std::vector<uint16_t> addBitsPos(positionPopulation.begin(), positionPopulation.begin() + numChangeBits);
		std::sort(addBitsPos.begin(), addBitsPos.end());
		std::vector<uint16_t>::iterator addBitsPosIt = addBitsPos.begin();

		IntChromosome mask = 0;
		IntChromosome tmp;
		int8_t totalShift = this->unusedBits - 1; // 0 based -- position 0 is the least significant bit
		int8_t trailingZeros = totalShift; // 1 based -- the value 1 means 1 trailing zero
		
		uint16_t onesCount = 0, zerosCount = 0;
		
		for(uint16_t i = 0; (i < this->numParts) && (addBitsPosIt != addBitsPos.end()); ++i) {
			tmp = this->chromosomeParts[i];
			do {
				tmp >>= trailingZeros + 1;
				
				trailingZeros = this->ctz(tmp); // get number of 0's before the first 1 (starting with least significant bit)
				
				if(trailingZeros >= Chromosome::BITS_PER_PART - totalShift - 1) {
					trailingZeros = Chromosome::BITS_PER_PART - totalShift - 1;
				}
				
				totalShift += trailingZeros + 1;
				
				// There are trailingZeros 0-bits between the last found 1 and this 1
				// Some of them may be toggled
				zerosCount += trailingZeros;
				while(addBitsPosIt != addBitsPos.end() && zerosCount > (*addBitsPosIt)) {
					mask |= (((IntChromosome) 1) << (((*addBitsPosIt) + onesCount + this->unusedBits) - i * Chromosome::BITS_PER_PART));
					++addBitsPosIt;
				}
				
				if(totalShift < Chromosome::BITS_PER_PART) {
					++onesCount;
				}
			} while(totalShift < Chromosome::BITS_PER_PART);
			
			this->chromosomeParts[i] ^= mask; // toggle all bits set in the mask
			
			totalShift = -1;
			trailingZeros = -1;
			mask = 0;
		}
	}
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << "After mutation: " << *this << std::endl;
	}
#endif
	
	return true;
}


std::ostream& operator<<(std::ostream &os, const Chromosome &ch) {
	ch.printBits(os, ch.chromosomeParts[0], ch.unusedBits);

	for (uint16_t i = 1; i < ch.numParts; ++i) {
		os << ' ';
		ch.printBits(os, ch.chromosomeParts[i]);
	}

	return os;
}

inline std::ostream& Chromosome::printBits(std::ostream &os, IntChromosome bits, uint16_t leaveOut) const {
	IntChromosome mask = ((IntChromosome) 1) << leaveOut;
	uint8_t delCount = 0;

	do {
		os << (((bits & mask) > 0) ? '1' : '0');
		if((++delCount % DELIMITER_POSITION) == 0) {
			delCount = 0;
			os << ' ';
		}
		mask <<= 1;
	} while(mask > 0);

	return os;
}

Rcpp::LogicalVector Chromosome::toLogicalVector() const {
	LogicalVector varVector;
	IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;

	for (uint16_t i = 0; i < this->numParts; ++i) {
		do {
			varVector.push_back((this->chromosomeParts[i] & mask) > 0);
			mask <<= 1;
		} while(mask > 0);
		mask = (IntChromosome) 1;
	}

	return varVector;
}

arma::uvec Chromosome::toColumnSubset() const {
	arma::uvec columnSubset(this->popcount());
	IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;
	uint16_t csIndex = 0;
	arma::uword truePos = 0;

	for (uint16_t i = 0; i < this->numParts && csIndex < columnSubset.n_elem; ++i) {
		do {
			if((this->chromosomeParts[i] & mask) > 0) {
				columnSubset[csIndex++] = truePos;
			}
			++truePos;
			mask <<= 1;
		} while(mask > 0 && csIndex < columnSubset.n_elem);
		mask = (IntChromosome) 1;
	}

	return columnSubset;

}

bool Chromosome::operator==(const Chromosome &ch) const {
	if(this->ctrl.chromosomeSize == ch.ctrl.chromosomeSize && this->numParts == ch.numParts) {
		for(uint16_t i = 0; i < this->numParts; ++i) {
			if(this->chromosomeParts[i] != ch.chromosomeParts[i]) {
				return false;
			}
		}
		return true;
	}
	return false;
}

bool Chromosome::operator!=(const Chromosome &ch) const {
	return !(*this == ch);
}

bool Chromosome::isFitterThan(const Chromosome &ch) const {
	if(this->fitness > ch.fitness) {
		return true;
	} else if(this->fitness < ch.fitness) {
		return false;
	}

	// Both chromosomes have (almost) same fitness
	// so check if this chormosome has less bits set than the other chromosome
	return (this->popcount() < ch.popcount());
}

Chromosome& Chromosome::operator=(const Chromosome &ch) {
	this->copyFrom(ch, true);
	return *this;
}

/**
 * Calculate the number of set bits in the chromosome
 * see the Wikipedia entry for "Hamming Weight"
 */
inline uint16_t Chromosome::popcount() const {
	uint16_t count = 0;

#ifdef HAVE_BUILTIN_POPCOUNT
	for(uint16_t i = 0; i < this->numParts; ++i) {
		count += __builtin_popcountl(this->chromosomeParts[i]);
	}
#else
	uint16_t i = 1;
	IntChromosome tmp = this->chromosomeParts[0] & ( INT_CHROMOSOME_MAX << this->unusedBits ); // "remove" first unused bits

	do {
		tmp -= (tmp >> 1) & Chromosome::M1;								// put count of each 2 bits into those 2 bits
		tmp = (tmp & Chromosome::M2) + ((tmp >> 2) & Chromosome::M2);	// put count of each 4 bits into those 4 bits
		tmp = (tmp + (tmp >> 4)) & Chromosome::M4;						// put count of each 8 bits into those 8 bits

		count += (tmp * Chromosome::H01) >> 56;							// adds left 8 bits of tmp + (tmp << 8) + (tmp << 16) + (tmp << 24) + ...

		tmp = this->chromosomeParts[i];
	} while(i++ < this->numParts);
#endif

	return count;
}

// returns BITS_PER_PART if no bit is set in mask
// positions start at 0
inline uint16_t Chromosome::ctz(IntChromosome mask) const {
// Only handle the case of 32 and 64 bit integers -- otherwise undefined behaviour
	if(mask == 0) {
		return Chromosome::BITS_PER_PART;
	}
	
#if defined HAVE_GCC_CTZ
	if(sizeof(IntChromosome) == sizeof(unsigned int)) {
		return  __builtin_ctz(mask);
	} else if(sizeof(IntChromosome) == sizeof(unsigned long)) {
		return  __builtin_ctzl(mask);
	} /*else if(sizeof(IntChromosome) == sizeof(unsigned long long)) {
		return  __builtin_ctzll(mask);
	}*/
	return Chromosome::BITS_PER_PART;
#elif defined HAVE_FFSL
	if(sizeof(IntChromosome) == sizeof(unsigned int)) {
		return  ffs(mask) - 1;
	} else if(sizeof(IntChromosome) == sizeof(unsigned long)) {
		return  ffsl(mask) - 1;
	}
	return Chromosome::BITS_PER_PART;
//#elif FFS_VS2005
//	unsigned long index;
//	unsigned char isNull = ((sizeof(mask) == 4) ? _BitScanForward(&index, mask) : _BitScanForward64(&index, mask));
//	return ((isNull == 0) ? (Chromosome::BITS_PER_PART) : (index - 1));
#else
	// Simple implementation of the "find first set" problem without loop	
	uint16_t index = 0;

	if (sizeof(mask) == 8) {
		if((mask & 0x00000000FFFFFFFF) == 0) {
			index += 32;
			mask >>= 32;
		}
	}

	if((mask & 0x0000FFFF) == 0) {
		index += 16;
		mask >>= 16;
	}
	if((mask & 0x000000FF) == 0) {
		index += 8;
		mask >>= 8;
	}
	if((mask & 0x0000000F) == 0) {
		index += 4;
		mask >>= 4;
	}
	if((mask & 0x00000003) == 0) {
		index += 2;
		mask >>= 2;
	}
	if((mask & 0x00000001) == 0) {
		index += 1;
	}
	
	return index;
#endif
}


inline IntChromosome Chromosome::runif() const {
	// Calculate the number of random numbers to combine for one IntChromosome
	static int8_t partsPerRand = (sizeof(IntChromosome) * BITS_PER_BYTE / RNG_MAX_BITS);
	
	if(partsPerRand > 1) {
		IntChromosome rand = 0;

		for(int8_t i = partsPerRand; i >= 0; --i) {
			rand |= (((IntChromosome) (INT_RNG_MAX * this->unifGen())) << (i * RNG_MAX_BITS));
		}

		return rand;
	} else {
		return (IntChromosome) (INT_CHROMOSOME_MAX * this->unifGen());
	}
}


inline void Chromosome::shuffle(std::vector<uint16_t>& pop, uint16_t fillLength, uint16_t shuffleLength) const {
	// Fill population with correct values
	if(shuffleLength == 1) {
		pop[0] = this->unifGen() * fillLength;
	} else {
		uint16_t i = 0;
		for(; i < fillLength; ++i) {
			pop[i] = i;
		}

		uint16_t randPos = 0;
		// Now shuffle population
		for(i = 0; i < shuffleLength; ++i) {
			randPos = i + this->unifGen() * (fillLength - i);
			std::swap(pop[i], pop[randPos]);
		}
	}
}

