#ifndef GenAlgPLS_Chromosome_h
#define GenAlgPLS_Chromosome_h

#include "config.h"

#include <inttypes.h>
#include <vector>
#include <exception>
#include <iostream>
#include <RcppArmadillo.h>

#include "Control.h"
#include "TruncatedGeomGenerator.h"
#include "VariablePositionPopulation.h"
#include "UnifGenerator__0__1.h"

class InvalidCopulationException : public Rcpp::exception {

public:
		InvalidCopulationException(const char *file, const int line);
};

class Chromosome {

public:
	Chromosome(const Control &ctrl, VariablePositionPopulation &varPosPop);
	Chromosome(const Chromosome &other, bool copyChromosomeParts = true);
//	~Chromosome();

	void mutate();
	std::vector<Chromosome> mateWith(const Chromosome &other);

	void setFitness(double fitness) { this->fitness = fitness; };
	double getFitness() const { return this->fitness; };

	Rcpp::LogicalVector toLogicalVector() const;
	arma::uvec toColumnSubset() const;

	bool isFitterThan(const Chromosome &ch) const;

	bool operator==(const Chromosome &ch) const;
	bool operator!=(const Chromosome &ch) const;
	Chromosome operator=(const Chromosome &ch) const;

	friend std::ostream& operator<<(std::ostream &os, const Chromosome &ch);
private:
	static const uint64_t M1 = 0x5555555555555555; // binary: 010101010101... (1 zero, 1 one)
	static const uint64_t M2 = 0x3333333333333333; // binary: 001100110011... (2 zeros, 2 ones)
	static const uint64_t M4 = 0x0f0f0f0f0f0f0f0f; // binary: 000011110000... (4 zeros, 4 ones)
	static const uint64_t H01 = 0x0101010101010101; // the sum of 256 to the power of 0,1,2,3...
	static const uint8_t BITS_PER_PART = sizeof(IntChromosome) * BITS_PER_BYTE;

	const Control &ctrl;
	const TruncatedGeomGenerator tgeom;
	const Rcpp::stats::UnifGenerator__0__1 unifGen;

	VariablePositionPopulation &varPosPop;

	void shuffle(std::vector<uint16_t>& pop, const uint16_t fillLength, const uint16_t shuffleLength) const;

	/*
	 * Init the internal used chromosome parts completely random taking
	 * the minimum and maximum number of set bits specified by the
	 * control object into account
	 */
	void initChromosomeParts();

	/*
	 * R's RNG only returns between 25 and 32 random bits
	 * so two random numbers must be "glued" together to form
	 * a 64bit random number
	 */
	IntChromosome runif() const;

	std::ostream& printBits(std::ostream &os, IntChromosome bits, uint16_t leaveOut = 0) const;

	uint16_t popcount() const;
	/*
	 * count trailing zeros
	 */
	uint16_t ctz(IntChromosome mask) const;

	uint16_t numParts;
	uint16_t unusedBits;

	std::vector<IntChromosome> chromosomeParts;
	double fitness;


};

#endif