#ifndef GenAlgPLS_Chromosome_h
#define GenAlgPLS_Chromosome_h

#include "config.h"

#include <inttypes.h>
#include <vector>
#include <exception>
#include <iostream>
#include <RcppArmadillo.h>
#include <Rcpp/stats/random/runif.h>

#include "Control.h"

class InvalidCopulationException : public Rcpp::exception {

public:
		InvalidCopulationException(const char *file, const int line);

};

class Chromosome {

public:
	Chromosome(const Control &ctrl);
	Chromosome(const Chromosome &other, bool copyChromosomeParts = true);
	~Chromosome();
	
	void mutate();
	std::vector<Chromosome*> copulateWith(const Chromosome &other);

	void setFitness(double fitness);
	double getFitness() const;

	Rcpp::LogicalVector toLogicalVector() const;
	arma::uvec toColumnSubset() const;
	
	bool isFitterThan(const Chromosome &ch) const;
	
	bool operator==(const Chromosome &ch) const;
	bool operator!=(const Chromosome &ch) const;
	
	friend std::ostream& operator<<(std::ostream &os, const Chromosome &ch);
private:
	static const uint_fast64_t M1 = 0x5555555555555555; // binary: 010101010101... (1 zero, 1 one)
	static const uint_fast64_t M2 = 0x3333333333333333; // binary: 001100110011... (2 zeros, 2 ones)
	static const uint_fast64_t M4 = 0x0f0f0f0f0f0f0f0f; // binary: 000011110000... (4 zeros, 4 ones)
	static const uint_fast64_t H01 = 0x0101010101010101; // the sum of 256 to the power of 0,1,2,3...

	const Rcpp::stats::UnifGenerator__0__1 unifGen;
	const uint16_t bitsPerPart;
	const Control ctrl;

	/*
	 * Generate a `size` numbers of at integral types with bits set at random.
	 * If onesRatio is less than 0.9, only that percentage of bits is set
	 * (actually, the number of set bits is binomial distributed with p = onesRatio)
	 */
	void generateRandomBits(uint_fast64_t* bits, uint16_t size, double onesRatio = 0.5) const;
		
	std::ostream& printBits(std::ostream &os, uint_fast64_t &bits, uint16_t leaveOut = 0) const;
	
	uint16_t popcount() const;
	
	uint16_t numParts;
	uint16_t unusedBits;
	
	uint_fast64_t* chromosomeParts;
	double fitness;
};

#endif