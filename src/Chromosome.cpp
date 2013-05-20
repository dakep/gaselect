//
//  Chromosome.cpp
//  
//

#include "config.h"

#include <inttypes.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <set>
#include <RcppArmadillo.h>
#include <Rcpp/stats/random/rbinom.h>

#include "Chromosome.h"
#include "GenAlg.h"

using namespace Rcpp;

InvalidCopulationException::InvalidCopulationException(const char *file, const int line) :
	Rcpp::exception("The two chromosomes are not compatible for mating", file, line) {
}

Chromosome::Chromosome(const Control &ctrl, VariablePositionPopulation &varPosPop) : ctrl(ctrl), tgeom(ctrl.mutationProbability), varPosPop(varPosPop) {
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
		
		if(part == 0)
			
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

std::vector<Chromosome> Chromosome::copulateWith(const Chromosome &other) {
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

void Chromosome::mutate() {
#ifdef TIMING_BENCHMARK
	timeval start, end;
	
	gettimeofday(&start, NULL);
#endif
	// The algorithm first randomly picks the number of bits that will be unset and the
	// number of bits that will be set according to a truncated geometric distribution
	// Then the positions for (un)setting the bits are drawn and finally these bits
	// are flipped.

	static std::vector<uint16_t> positionPopulation(this->ctrl.chromosomeSize);
	
	uint16_t currentlySetBits = this->popcount();
	uint16_t currentlyUnsetBits = this->ctrl.chromosomeSize - currentlySetBits;

	uint16_t numAddBits = 0, numRemoveBits = 0;
	
	if(this->ctrl.minVariables - currentlySetBits > 0) { // Too few bits are set -- we MUST add some
		numAddBits = (this->ctrl.minVariables - currentlySetBits);
	} else if(currentlySetBits - this->ctrl.maxVariables > 0) {
		numRemoveBits = currentlySetBits - this->ctrl.maxVariables;
	}
	
	
	if(this->ctrl.maxVariables - currentlySetBits > 0) { // We may add some bits
		numAddBits += this->tgeom(this->ctrl.maxVariables - currentlySetBits - numAddBits);
	}

	if(currentlySetBits - this->ctrl.minVariables > 0) { // We may remove some bits
		numRemoveBits += this->tgeom(currentlySetBits - this->ctrl.minVariables - numRemoveBits);
	}

#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << "Adding " << numAddBits << " variables and removing " << numRemoveBits << " variables" << std::endl;
	}
#endif

	if(numAddBits == 0 && numRemoveBits == 0) {
		return;
	}
	
	uint16_t randPos = 0;
	
	// Choose the fastest way to mutate the chromosome based
	// on the number of set/unset bits and the number of bits
	// to add/remove
	
	if(numRemoveBits == 0) {
		// Bits only have to be added

		if(numAddBits == 1 && currentlyUnsetBits * RATIO_RANDOM_SEARCH > currentlySetBits) {
			// If the number of bits to be added is 1 and the majority of the bits is 0, a random search may be faster than the method below
			uint16_t part = 0;
			IntChromosome mask = 0;
			do {
				randPos = this->unusedBits + this->unifGen() * this->ctrl.chromosomeSize;
				part = randPos / Chromosome::BITS_PER_PART;
				mask = (((IntChromosome) 1) << (randPos % Chromosome::BITS_PER_PART));
			} while ((this->chromosomeParts[part] & mask) > 0); // Exit if the bit at the random position is not yet 1
		
			this->chromosomeParts[part] |= mask;
		} else {
			// Same algorithm as for general case but removed all code for removing bits
			this->shuffle(positionPopulation, currentlyUnsetBits, numAddBits);
			std::vector<uint16_t> addBitsPos(positionPopulation.begin(), positionPopulation.begin() + numAddBits);			
			std::sort(addBitsPos.begin(), addBitsPos.end());
			std::vector<uint16_t>::iterator addBitsPosIt = addBitsPos.begin();
			IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;
			uint16_t zerosCount = 0;
			
			for(uint16_t i = 0; i < this->numParts && addBitsPosIt != addBitsPos.end(); ++i) {
				do {
					if((this->chromosomeParts[i] & mask) == 0) { // bit is 0
						if(addBitsPosIt != addBitsPos.end() && zerosCount++ == (*addBitsPosIt)) {
							this->chromosomeParts[i] ^= mask;
							++addBitsPosIt;
						}
					}
					
					mask <<= 1;
				} while (mask > 0 && addBitsPosIt != addBitsPos.end());
				mask = (IntChromosome) 1;
			}
			
		}
	} else if(numAddBits == 0) {
		// Bits only have to be removed
		if(numRemoveBits == 1 && currentlySetBits * RATIO_RANDOM_SEARCH > currentlyUnsetBits) {
			// If the number of bits to be removed is small and the majority of the bits is 1, a random search may be faster than the method below
			uint16_t part = 0;
			IntChromosome mask = 0;
			do {
				randPos = this->unusedBits + this->unifGen() * this->ctrl.chromosomeSize;
				part = randPos / Chromosome::BITS_PER_PART;
				mask = (((IntChromosome) 1) << (randPos % Chromosome::BITS_PER_PART));
			} while ((this->chromosomeParts[part] & mask) == 0); // Exit if the bit at the random position is not yet 0
			
			this->chromosomeParts[part] &= ~mask; // Remove only the bit at the random position
			
		} else {
			this->shuffle(positionPopulation, currentlySetBits, numRemoveBits);
			std::vector<uint16_t> removeBitsPos(positionPopulation.begin(), positionPopulation.begin() + numRemoveBits);
			std::sort(removeBitsPos.begin(), removeBitsPos.end());
			std::vector<uint16_t>::iterator removeBitsPosIt = removeBitsPos.begin();
			IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;
			uint16_t onesCount = 0;
			
			for(uint16_t i = 0; i < this->numParts && removeBitsPosIt != removeBitsPos.end(); ++i) {
				do {
					if((this->chromosomeParts[i] & mask) > 0) { // bit is 1
						if(removeBitsPosIt != removeBitsPos.end() && onesCount++ == (*removeBitsPosIt)) {
							this->chromosomeParts[i] ^= mask;
							++removeBitsPosIt;
						}
					}
					
					mask <<= 1;
				} while (mask > 0 && removeBitsPosIt != removeBitsPos.end());
				mask = (IntChromosome) 1;
			}
		}
	} else {
		// General case
		
		this->shuffle(positionPopulation, currentlyUnsetBits, numAddBits);
		std::vector<uint16_t> addBitsPos(positionPopulation.begin(), positionPopulation.begin() + numAddBits);

		this->shuffle(positionPopulation, currentlySetBits, numRemoveBits);
		std::vector<uint16_t> removeBitsPos(positionPopulation.begin(), positionPopulation.begin() + numRemoveBits);
		
		std::sort(addBitsPos.begin(), addBitsPos.end());
		std::sort(removeBitsPos.begin(), removeBitsPos.end());
		
		std::vector<uint16_t>::iterator addBitsPosIt = addBitsPos.begin();
		std::vector<uint16_t>::iterator removeBitsPosIt = removeBitsPos.begin();
		
		IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;
		
		uint16_t onesCount = 0, zerosCount = 0;
		
		for(uint16_t i = 0; i < this->numParts && (addBitsPosIt != addBitsPos.end() || removeBitsPosIt != removeBitsPos.end()); ++i) {
			do {
				if((this->chromosomeParts[i] & mask) > 0) { // bit is 1
					if(removeBitsPosIt != removeBitsPos.end() && onesCount++ == (*removeBitsPosIt)) {
						this->chromosomeParts[i] ^= mask;
						++removeBitsPosIt;
					}
				} else { // bit is 0
					if(addBitsPosIt != addBitsPos.end() && zerosCount++ == (*addBitsPosIt)) {
						this->chromosomeParts[i] ^= mask;
						++addBitsPosIt;
					}
				}
				
				mask <<= 1;
			} while (mask > 0 && (addBitsPosIt != addBitsPos.end() || removeBitsPosIt != removeBitsPos.end()));
			mask = (IntChromosome) 1;
		}
	}
#ifdef ENABLE_DEBUG_VERBOSITY
	if(this->ctrl.verbosity == DEBUG_VERBOSE) {
		Rcout << "After mutation: " << *this << std::endl;
	}
#endif

#ifdef TIMING_BENCHMARK

	gettimeofday(&end, NULL);
	Rcout << "Mutation took " << (end.tv_sec * 1000.0 + (end.tv_usec / 1000.0)) - (start.tv_sec * 1000.0 + (start.tv_usec / 1000.0)) << " milliseconds" << std::endl;

#endif
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

Chromosome Chromosome::operator=(const Chromosome &ch) const {
	return Chromosome(ch);
}

/**
 * Calculate the number of set bits in the chromosome
 * see the Wikipedia entry for "Hamming Weight"
 */
inline uint16_t Chromosome::popcount() const {
	uint16_t count = 0;
	
#ifdef USE_BUILTIN_POPCOUNT
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


inline IntChromosome Chromosome::runif() const {
#if RANDS_PER_INT_CHROMOSOME > 1
	IntChromosome rand = 0;
	
	for(uint16_t i = 0; i < RANDS_PER_INT_CHROMOSOME; ++i) {
		rand |= (((IntChromosome) (INT_RNG_MAX * this->unifGen())) << (i * RNG_MAX_BITS));
	}

	return rand;
#else
	return (IntChromosome) (INT_CHROMOSOME_MAX * this->unifGen());
#endif
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

