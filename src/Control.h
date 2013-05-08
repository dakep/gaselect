//
//  Control.h
//  
//

#ifndef GenAlgPLS_Control_h
#define GenAlgPLS_Control_h

#include "config.h"

#include <iostream>
#include <inttypes.h>
#include <RcppArmadillo.h>

enum VerbosityLevel {
	OFF = 0,
	ON,
	MORE_VERBOSE,
	DEBUG_VERBOSE
};

class Control {

public:
	Control() : chromosomeSize(0), popSize(0), numGenerations(0), elitism(0), onesRatio(0.0), mutate0To1Probability(0.0), mutate1To0Probability(0.0), verbosity(OFF) {}
	
	void setChromosomeSize(uint16_t chromosomeSize) { this->chromosomeSize = chromosomeSize; };
	uint16_t getChromosomeSize() const { return this->chromosomeSize; };

	void setPopulationSize(uint16_t popSize);
	
	uint16_t getPopulationSize() const { return this->popSize; };
	
	void setNumberOfGenerations(uint16_t numGenerations) { this->numGenerations = numGenerations; };
	uint16_t getNumberOfGenerations() const { return this->numGenerations; };
	
	void setElitism(uint16_t elitism) { this->elitism = elitism; };
	uint16_t getElitism() const { return this->elitism; }
	
	void setOnesRatio(double onesRatio);
	double getOnesRatio() const { return this->onesRatio; }
	
	void setMutate0To1Probability(double mutate0To1Probability);
	double getMutate0To1Probability() const { return this->mutate0To1Probability; };
	double getMutate1To0Probability() const { return this->mutate1To0Probability; };
	
	void setVerbosity(enum VerbosityLevel verbosity) { this->verbosity = verbosity; };
	enum VerbosityLevel getVerbosity() const { return this->verbosity; };
	
	friend std::ostream& operator<<(std::ostream &os, const Control &ctrl);
private:
	void setMutate1To0Probability();
	
	uint16_t chromosomeSize;
	uint16_t popSize;
	uint16_t numGenerations;
	uint16_t elitism;
	double onesRatio;
	
	double mutate0To1Probability;
	double mutate1To0Probability;
	
	enum VerbosityLevel verbosity;
};

#endif
