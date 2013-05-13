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

Chromosome::Chromosome(const Control &ctrl, VariablePositionPopulation &varPosPop) : ctrl(ctrl), varPosPop(varPosPop) {
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
	this->chromosomeParts = new IntChromosome[this->numParts];
	this->initChromosomeParts();
}

Chromosome::Chromosome(const Chromosome &other, bool copyChromosomeParts) : ctrl(other.ctrl), varPosPop(other.varPosPop) {
	this->fitness = other.fitness;
	this->numParts = other.numParts;
	this->unusedBits = other.unusedBits;
	
	this->chromosomeParts = new IntChromosome[this->numParts];
	
	// Copy chromosome parts
	if(copyChromosomeParts) {
		std::copy(other.chromosomeParts, other.chromosomeParts + this->numParts, this->chromosomeParts);
	}	
}

Chromosome::~Chromosome() {
	delete [] this->chromosomeParts;
}

void Chromosome::setFitness(double fitness) {
	this->fitness = fitness;
}

double Chromosome::getFitness() const {
	return this->fitness;
}

void Chromosome::mutate() {
#ifdef TIMING_BENCHMARK
	timeval start, end;
	
	gettimeofday(&start, NULL);
#endif

	// The algorithm first randomly picks the number of 1's resp. 0's that are to be flipped
	// The next step is to randomly pick the position of these flips
	// In the last step, those positions are flipped
		
	uint16_t currentlySetBits = this->popcount();
	uint16_t currentlyUnsetBits = this->ctrl.chromosomeSize - currentlySetBits;
	
	uint16_t flip1To0Count = this->unifGen() * (this->ctrl.maxVariables - currentlySetBits);
	uint16_t flip0To1Count = this->unifGen() * (this->ctrl.maxVariables - currentlySetBits + flip1To0Count);
	
	// Unset numBitsToUnset
	
	// Set flip1To0Count

	std::vector<uint16_t> flip1To0Pos(flip1To0Count);
	std::vector<uint16_t> flip0To1Pos(flip0To1Count);
		
	uint16_t randPos = 0;
	uint16_t i = 0;
	
	std::vector<uint16_t> posPop;
	posPop.reserve((currentlyUnsetBits < currentlySetBits) ? currentlySetBits : currentlyUnsetBits);
	
	// The new position is somewhere between i and # of 1's,
	// so the random number can be between 0 and # of 1's - i;
	uint16_t maxPos = currentlySetBits;

	for(; i < currentlySetBits; ++i) {
		posPop[i] = i;
	}

	for(i = 0; i < flip1To0Count; ++i) {
		randPos = i + (uint16_t) maxPos-- * this->unifGen();
		flip1To0Pos[i] = posPop[randPos];
		std::swap(posPop[randPos], posPop[i]);
	}
	
	for(; i < currentlyUnsetBits; ++i) {
		posPop[i] = i;
	}

	maxPos = currentlyUnsetBits;
	
	for(i = 0; i < flip0To1Count; ++i) {
		randPos = i + (uint16_t) maxPos-- * this->unifGen();
		flip0To1Pos[i] = posPop[randPos];
		std::swap(posPop[randPos], posPop[i]);
	}
	
	std::sort(flip1To0Pos.begin(), flip1To0Pos.end());
	std::sort(flip0To1Pos.begin(), flip0To1Pos.end());

	std::vector<uint16_t>::iterator flip1To0PosIt = flip1To0Pos.begin();
	std::vector<uint16_t>::iterator flip0To1PosIt = flip0To1Pos.begin();
	
	IntChromosome mask = ((IntChromosome) 1) << this->unusedBits;
	
	uint16_t onesCount = 0, zerosCount = 0;
	
	for(i = 0; i < this->numParts; ++i) {
		do {
			if(flip1To0Count > 0 && (this->chromosomeParts[i] & mask) > 0) { // bit is 1
				if(onesCount++ == (*flip1To0PosIt)) {
					this->chromosomeParts[i] ^= mask;
					++flip1To0PosIt;
				}
			} else { // bit is 0
				if(flip0To1Count > 0 && zerosCount++ == (*flip0To1PosIt)) {
					this->chromosomeParts[i] ^= mask;
					++flip0To1PosIt;
				}
			}
			
			mask <<= 1;
		} while (mask > 0);
		mask = (IntChromosome) 1;
	}
	
// Slower, but simpler algorithm
	
//	uint16_t i = 0;
//	int j = Chromosome::BITS_PER_PART - this->unusedBits - 1;
//	IntChromosome checkMask = ((IntChromosome) 1) << j;
//	double prob1To0 = ((1 - this->ctrl.getOnesRatio()) / this->ctrl.getOnesRatio()) * this->ctrl.getMutationProbability();
//	double rand = 0.0;
//	
//	for(; i < this->numParts; ++i) {
//		for(; j >= 0; --j) {
//			//			varVector.push_back((this->chromosomeParts[i] & mask) > 0);
//			
//			rand = this->unifGen();
//			
//			if((this->chromosomeParts[i] & checkMask) > 0) { // Bit is 1
//				if(this->unifGen() < prob1To0) {
//					this->chromosomeParts[i] ^= checkMask;
//				}
//			} else { // Bit is zero
//				if(this->unifGen() < this->ctrl.getMutationProbability()) { // Should it be flipped?
//					this->chromosomeParts[i] ^= checkMask;
//				}
//			}
//			
//			checkMask >>= 1;
//		}
//		j = Chromosome::BITS_PER_PART - 1;
//	}

#ifdef TIMING_BENCHMARK

	gettimeofday(&end, NULL);
	Rcout << "Mutation took " << (end.tv_sec * 1000.0 + (end.tv_usec / 1000.0)) - (start.tv_sec * 1000.0 + (start.tv_usec / 1000.0)) << " milliseconds" << std::endl;

#endif
}

inline void Chromosome::initChromosomeParts() {
	uint16_t bitsToSet = this->ctrl.minVariables + this->unifGen() * (this->ctrl.minVariables - this->ctrl.maxVariables);
	
	VariablePositionPopulation::const_iterator setPosIter = this->varPosPop.shuffle(bitsToSet);

	uint16_t part = 0;
	uint16_t offset = 0;
	
	for(; setPosIter != this->varPosPop.end(); ++setPosIter) {
		part = (*setPosIter) / Chromosome::BITS_PER_PART;
		offset = (*setPosIter) % Chromosome::BITS_PER_PART;
		
		this->chromosomeParts[part] |= (((IntChromosome) 1) << offset);
	}		
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
			if(this->ctrl.verbosity == MORE_VERBOSE) {
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
				this->printBits(Rcout, child1->chromosomeParts[i], (i == 0) ? this->unusedBits : 0) << std::endl
				<< "Child 2:";
				this->printBits(Rcout, child2->chromosomeParts[i], (i == 0) ? this->unusedBits : 0) << std::endl;
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

std::ostream& operator<<(std::ostream &os, const Chromosome &ch) {
	ch.printBits(os, ch.chromosomeParts[0], ch.unusedBits);
	
	for (uint16_t i = 1; i < ch.numParts; ++i) {
		os << ' ';
		ch.printBits(os, ch.chromosomeParts[i]);
	}
	
	return os;
}

inline std::ostream& Chromosome::printBits(std::ostream &os, IntChromosome &bits, uint16_t leaveOut) const {
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
		rand |= (((IntChromosome) (INT_CHROMOSOME_MAX * this->unifGen())) << (i * RNG_MAX_BITS));
	}
	
	return rand;
#else
	return (IntChromosome) (INT_CHROMOSOME_MAX * this->unifGen());
#endif
}


