//
//  Control.cpp
//  
//

#include "Control.h"

void Control::setPopulationSize(uint16_t popSize) {
	if((popSize % 2) > 0) {
		this->popSize = popSize - 1;
	} else {
		this->popSize = popSize;
	}
}

void Control::setOnesRatio(double onesRatio) {
	this->onesRatio = onesRatio;

	if(this->mutate0To1Probability > 0) {
		this->setMutate1To0Probability();
	}
}

void Control::setMutate0To1Probability(double mutate0To1Probability) {
	this->mutate0To1Probability = mutate0To1Probability;
	
	if(this->onesRatio > 0) {
		this->setMutate1To0Probability();
	}
}

inline void Control::setMutate1To0Probability() {
	this->mutate1To0Probability = ((1.0 - this->onesRatio) / this->onesRatio) * this->mutate0To1Probability;
}

std::ostream& operator<<(std::ostream &os, const Control &ctrl) {
	os << "Chromosome size: " << ctrl.chromosomeSize << std::endl
		<< "Population size: " << ctrl.popSize << std::endl
		<< "Number of generations: " << ctrl.numGenerations << std::endl
		<< "Number of elite chromosomes to keep: " << ctrl.elitism << std::endl
		<< "1 to 0 ratio: " << ctrl.onesRatio << std::endl
		<< "Probability to mutate 0 to 1: " << ctrl.mutate0To1Probability << std::endl
	<< "Verbosity Level: " << ctrl.verbosity << std::endl;

	return os;
}
