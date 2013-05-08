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

Chromosome::Chromosome(const Control &ctrl) : unifGen(), bitsPerPart(sizeof(uint_fast64_t) * BITS_PER_BYTE), ctrl(ctrl) {
	// Determine the number of uint_fast64_t bit values that are
	// needed to represent all genes
	this->numParts = (uint16_t) this->ctrl.getChromosomeSize() / this->bitsPerPart;

	int rest = this->ctrl.getChromosomeSize() % this->bitsPerPart;
	
	if(rest > 0) {
		this->numParts += 1;
		this->unusedBits = this->bitsPerPart - rest;
	} else {
		this->unusedBits = 0;
	}

	this->fitness = 0.0;

	// Initialize chromosome parts randomly
	this->chromosomeParts = new uint_fast64_t[this->numParts];
	this->generateRandomBits(this->chromosomeParts, this->numParts, this->ctrl.getOnesRatio());
}

Chromosome::Chromosome(const Chromosome &other, bool copyChromosomeParts) : bitsPerPart(sizeof(uint_fast64_t) * BITS_PER_BYTE), ctrl(other.ctrl) {
	this->fitness = other.fitness;
	this->numParts = other.numParts;
	this->unusedBits = other.unusedBits;
	
	this->chromosomeParts = new uint_fast64_t[this->numParts];
	
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
		
	uint16_t numSetBits = this->popcount();
	uint16_t numNotSetBits = this->ctrl.getChromosomeSize() - numSetBits;
	
	// numFlip1To0Gen ~ Bin(# of 1's; ((1 - p) / p) * mutProb)
	Rcpp::stats::BinomGenerator numFlip1To0Gen(numSetBits, this->ctrl.getMutate1To0Probability());
	// numFlip0To1Gen ~ Bin(# of 0's, mutProb)
	Rcpp::stats::BinomGenerator numFlip0To1Gen(numNotSetBits, this->ctrl.getMutate0To1Probability());
	
	uint16_t flip1To0Count = numFlip1To0Gen();
	uint16_t flip0To1Count = numFlip0To1Gen();

	std::vector<uint16_t> flip1To0Pos(flip1To0Count);
	std::vector<uint16_t> flip0To1Pos(flip0To1Count);
		
	uint16_t randPos = 0;
	uint16_t i = 0;
		
	uint16_t *posPop = new uint16_t[(numNotSetBits < numSetBits) ? numSetBits : numNotSetBits];
	
	// The new position is somewhere between i and # of 1's,
	// so the random number can be between 0 and # of 1's - i;
	uint16_t maxPos = numSetBits;

	for(; i < numSetBits; ++i) {
		posPop[i] = i;
	}

	for(i = 0; i < flip1To0Count; ++i) {
		randPos = i + (uint16_t) maxPos-- * (this->unifGen() - DISCRETE_CORRECTION);
		flip1To0Pos[i] = posPop[randPos];
		std::swap(posPop[randPos], posPop[i]);
	}
	
	for(; i < numNotSetBits; ++i) {
		posPop[i] = i;
	}

	maxPos = numNotSetBits;
	
	for(i = 0; i < flip0To1Count; ++i) {
		randPos = i + (uint16_t) maxPos-- * (this->unifGen() - DISCRETE_CORRECTION);
		flip0To1Pos[i] = posPop[randPos];
		std::swap(posPop[randPos], posPop[i]);
	}
	
	std::sort(flip1To0Pos.begin(), flip1To0Pos.end());
	std::sort(flip0To1Pos.begin(), flip0To1Pos.end());

	std::vector<uint16_t>::iterator flip1To0PosIt = flip1To0Pos.begin();
	std::vector<uint16_t>::iterator flip0To1PosIt = flip0To1Pos.begin();
	
	uint_fast64_t mask = ((uint_fast64_t) 1) << this->unusedBits;
	
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
		mask = (uint_fast64_t) 1;
	}
	
// Slower, but simpler algorithm
	
//	uint16_t i = 0;
//	int j = this->bitsPerPart - this->unusedBits - 1;
//	uint_fast64_t checkMask = ((uint_fast64_t) 1) << j;
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
//		j = this->bitsPerPart - 1;
//	}

#ifdef TIMING_BENCHMARK

	gettimeofday(&end, NULL);
	Rcout << "Mutation took " << (end.tv_sec * 1000.0 + (end.tv_usec / 1000.0)) - (start.tv_sec * 1000.0 + (start.tv_usec / 1000.0)) << " milliseconds" << std::endl;

#endif
}

/*
 * Generate a stream of random bits with a specified probability of setting a bit to 1
 * This algorithm only needs as many iterations as bits are set (E[it] = chromosomeSize * onesRatio)
 * args:
 *		bits ... Array of x-bit sized integral values
 *		size ... Size of the array
 *		onesRatio ... Probability of setting a bit to 1
 */
void Chromosome::generateRandomBits(uint_fast64_t* bits, uint16_t size, double onesRatio) const {
	Rcpp::stats::BinomGenerator binomGen(this->bitsPerPart, onesRatio);

	uint_fast64_t randNr = 0;
	int bitsToSet = 0;
	bool negateResult = false;
	uint_fast64_t setBitAt = 0; // a number with only the bit set at a random position between 0 and e.g. 63
		
	for (uint16_t i = 0; i < size; ++i) {
		randNr = 0;
		bitsToSet = (int) binomGen();
		negateResult = (bitsToSet > (this->bitsPerPart / 2));
		setBitAt = 0;
		
#ifdef ENABLE_DEBUG_VERBOSITY
		if(this->ctrl.getVerbosity() == DEBUG_VERBOSE) {
			Rcout << "Setting " << bitsToSet << "/" << this->bitsPerPart << " bits" << std::endl;
		}
#endif
		
		if(negateResult) {
			// More bits have to be set than not set
			// so the bits that should be unset are set
			// and the bits will be negated in the end
			bitsToSet = this->bitsPerPart - bitsToSet;
		}
		
		for(; bitsToSet > 0; --bitsToSet) {
			
			// Search for a position that is not yet set
			// This may lead to a very long loop!
			do {
				setBitAt = ((uint_fast64_t) 1) << (uint_fast64_t)(this->bitsPerPart * (this->unifGen() - DISCRETE_CORRECTION));
			} while ((randNr & setBitAt) > 0);
			
			randNr |= setBitAt;
		}
		
		if (negateResult) {
			 bits[i] = ~randNr;
		} else {
			bits[i] = randNr;
		}
	}
}

std::vector<Chromosome*> Chromosome::copulateWith(const Chromosome &other) {
#ifdef TIMING_BENCHMARK
	timeval start, end;
	
	gettimeofday(&start, NULL);
#endif
	if(other.ctrl.getChromosomeSize() != this->ctrl.getChromosomeSize()) {
		throw InvalidCopulationException(__FILE__, __LINE__);
	}
	
	std::vector<Chromosome*> children;
	
	Chromosome* child1 = new Chromosome(*this, false);
	Chromosome* child2 = new Chromosome(*this, false);

	uint_fast64_t randomMask = 0;
	uint_fast64_t negRandomMask = 0;
	
//	uint_fast64_t *randomMask = new uint_fast64_t[this->numParts];
//	this->generateRandomBits(randomMask, this->numParts);
	
	for(uint16_t i = 0; i < this->numParts; ++i) {
		if(this->chromosomeParts[i] == other.chromosomeParts[i]) {
			// Just copy the chromosome part to both children if it is the same
			// for both parents
			child1->chromosomeParts[i] = child2->chromosomeParts[i] = this->chromosomeParts[i];
			if(this->ctrl.getVerbosity() == MORE_VERBOSE) {
				Rcout << "Chromosome part is the same for both parents -- copy part to both children" << std::endl;
			}
		} else {
			// Randomly pick some bits from one chromosome and some bits from the other
			// chromosome
			do {
				this->generateRandomBits(&randomMask, 1);
				if(i == 0) {
					randomMask <<= this->unusedBits;
				}
			} while(randomMask == 0);
			negRandomMask = ~randomMask;
			
			child1->chromosomeParts[i] = (this->chromosomeParts[i] & randomMask) | (other.chromosomeParts[i] & negRandomMask);
			child2->chromosomeParts[i] = (this->chromosomeParts[i] & negRandomMask) | (other.chromosomeParts[i] & randomMask);

#ifdef ENABLE_DEBUG_VERBOSITY
			if(this->ctrl.getVerbosity() == DEBUG_VERBOSE) {
				Rcout << "Mask for part " << i << ": ";
				this->printBits(Rcout, randomMask, (i == 0) ? this->unusedBits : 0) << std::endl;
				
				Rcout << "Resulting parts:" << std::endl
				<< "Child 1:";
				this->printBits(Rcout, child1->chromosomeParts[i], (i == 0) ? this->unusedBits : 0) << std::endl
				<< "Child 2:";
				this->printBits(Rcout, child2->chromosomeParts[i], (i == 0) ? this->unusedBits : 0) << std::endl;
			}	
		}
	}
#endif
	
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

inline std::ostream& Chromosome::printBits(std::ostream &os, uint_fast64_t &bits, uint16_t leaveOut) const {
	uint_fast64_t mask = ((uint_fast64_t) 1) << leaveOut;
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
	uint_fast64_t mask = ((uint_fast64_t) 1) << this->unusedBits;
	
	for (uint16_t i = 0; i < this->numParts; ++i) {
		do {
			varVector.push_back((this->chromosomeParts[i] & mask) > 0);
			mask <<= 1;
		} while(mask > 0);
		mask = (uint_fast64_t) 1;
	}
	
	return varVector;
}

arma::uvec Chromosome::toColumnSubset() const {
	arma::uvec columnSubset(this->popcount());
	uint_fast64_t mask = ((uint_fast64_t) 1) << this->unusedBits;
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
		mask = (uint_fast64_t) 1;
	}
	
	return columnSubset;

}

bool Chromosome::operator==(const Chromosome &ch) const {
	if(this->ctrl.getChromosomeSize() == ch.ctrl.getChromosomeSize() && this->numParts == ch.numParts && this->bitsPerPart == ch.bitsPerPart) {
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

/**
 * Calculate the number of set bits in the chromosome
 * see the Wikipedia entry for "Hamming Weight"
 */
inline uint16_t Chromosome::popcount() const {
	uint16_t count = 0;
	uint16_t i = 0;
	uint_fast64_t tmp = this->chromosomeParts[0] & ( UINT_FAST64_MAX << this->unusedBits ); // "remove" first unused bits

	do {
		tmp -= (tmp >> 1) & Chromosome::M1;								// put count of each 2 bits into those 2 bits
		tmp = (tmp & Chromosome::M2) + ((tmp >> 2) & Chromosome::M2);	// put count of each 4 bits into those 4 bits
		tmp = (tmp + (tmp >> 4)) & Chromosome::M4;						// put count of each 8 bits into those 8 bits
		
		count += (tmp * Chromosome::H01) >> 56;							// adds left 8 bits of tmp + (tmp << 8) + (tmp << 16) + (tmp << 24) + ...
		
		tmp = this->chromosomeParts[i];
	} while(++i < this->numParts);
	
	return count;
}

