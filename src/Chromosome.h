#ifndef GenAlgPLS_Chromosome_h
#define GenAlgPLS_Chromosome_h

#include "config.h"

#ifdef HAVE_CLIMITS
#include <climits>
#elif HAVE_LIMITS_H
#include <limits.h>
#endif

#include <vector>
#include <exception>
#include <iostream>
#include <RcppArmadillo.h>

#include "Control.h"
#include "TruncatedGeomGenerator.h"
#include "VariablePositionPopulation.h"
#include "UnifGenerator_0_1.h"

class InvalidCopulationException : public Rcpp::exception {

public:
	InvalidCopulationException(const char *file, const int line) : Rcpp::exception("The two chromosomes are not compatible for mating", file, line) {};
};

class Chromosome {

public:
	Chromosome(const Control &ctrl, VariablePositionPopulation &varPosPop, UnifGenerator_0_1& unifGen, bool randomInit = true);
	Chromosome(const Chromosome &other, bool copyChromosomeParts = true);
//	~Chromosome();
	
	/**
	 * @return bool Returns true if mutation occurred, false otherwise
	 */
	bool mutate(UnifGenerator_0_1& unifGen);
	void mateWith(const Chromosome &other, UnifGenerator_0_1& unifGen, Chromosome& child1, Chromosome& child2);

	void setFitness(double fitness) { this->fitness = fitness; };
	double getFitness() const { return this->fitness; };

	Rcpp::LogicalVector toLogicalVector() const;
	arma::uvec toColumnSubset() const;

	bool isFitterThan(const Chromosome &ch) const;

	bool operator==(const Chromosome &ch) const;
	bool operator!=(const Chromosome &ch) const;
	Chromosome& operator=(const Chromosome &ch);
	
	uint16_t getVariableCount() const;

	friend std::ostream& operator<<(std::ostream &os, const Chromosome &ch);
private:
#if !(defined HAVE_BUILTIN_POPCOUNTLL | defined HAVE_BUILTIN_POPCOUNTL)
	static const IntChromosome M1 = 0x5555555555555555; // binary: 010101010101... (1 zero, 1 one)
	static const IntChromosome M2 = 0x3333333333333333; // binary: 001100110011... (2 zeros, 2 ones)
	static const IntChromosome M4 = 0x0f0f0f0f0f0f0f0f; // binary: 000011110000... (4 zeros, 4 ones)
	static const IntChromosome H01 = 0x0101010101010101; // the sum of 256 to the power of 0,1,2,3...
#endif
	static const uint8_t BITS_PER_PART = sizeof(IntChromosome) * BITS_PER_BYTE;

#ifdef INT_CHROMOSOME_MAX_VAL
	static const IntChromosome INT_CHROMOSOME_MAX = INT_CHROMOSOME_MAX_VAL;
#else
	static IntChromosome INT_CHROMOSOME_MAX;
	static IntChromosome getIntChromosomeMax();
#endif
	const Control &ctrl;
	const TruncatedGeomGenerator rtgeom;

	uint16_t numParts;
	uint16_t unusedBits;
	
	// Array with the chromosome parts
	// If not all bits are used, the k least significant bits of the 1st(!) part are not used (k = this->unusedBits)
	std::vector<IntChromosome> chromosomeParts;
	double fitness;

	void shuffle(std::vector<uint16_t>& pop, const uint16_t fillLength, const uint16_t shuffleLength, UnifGenerator_0_1& unifGen) const;

	/*
	 * Init the internal used chromosome parts completely random taking
	 * the minimum and maximum number of set bits specified by the
	 * control object into account
	 */
	void initChromosomeParts(UnifGenerator_0_1& unifGen, VariablePositionPopulation &varPosPop);

	/*
	 * R's RNG only returns between 25 and 32 random bits
	 * so two random numbers must be "glued" together to form
	 * a 64bit random number
	 */
	IntChromosome runif(UnifGenerator_0_1& unifGen) const;

	std::ostream& printBits(std::ostream &os, IntChromosome bits, uint16_t leaveOut = 0) const;

	/*
	 * count trailing zeros
	 */
	uint16_t ctz(IntChromosome mask) const;

	void copyFrom(const Chromosome& ch, bool copyChromosomeParts);

#if !(defined HAVE_BUILTIN_POPCOUNTLL | defined HAVE_BUILTIN_POPCOUNTL)
	uint16_t popcount(IntChromosome x) const {
		if(x == 0) {
			return 0;
		}
		
		x -= (x >> 1) & Chromosome::M1;								// put count of each 2 bits into those 2 bits
		x = (x & Chromosome::M2) + ((x >> 2) & Chromosome::M2);	// put count of each 4 bits into those 4 bits
		x = (x + (x >> 4)) & Chromosome::M4;						// put count of each 8 bits into those 8 bits
		return ((x * Chromosome::H01) >> 56);							// adds left 8 bits of tmp + (tmp << 8) + (tmp << 16) + (tmp << 24) + ...
	}
#endif


};

#endif