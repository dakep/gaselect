//
//  Control.cpp
//
//

#include "config.h"

#include <algorithm>
#include "Control.h"

Control::Control(const uint16_t chromosomeSize, const uint16_t popSize, const uint16_t numGenerations, const uint16_t elitism, const uint16_t minVariables, const uint16_t maxVariables, const uint16_t maxMatingTries, const double mutationProbability, const uint16_t numThreads, const enum VerbosityLevel verbosity) :
	chromosomeSize(chromosomeSize), populationSize(popSize), numGenerations(numGenerations), elitism(elitism), minVariables(minVariables), maxVariables(maxVariables),
	maxMatingTries(maxMatingTries), mutationProbability(mutationProbability), numThreads(numThreads), verbosity(verbosity) {

		this->variablePositionPopulation.reserve(this->chromosomeSize);

		uint16_t i = 0, j = 1;
		for(; j < this->chromosomeSize; i += 2, j += 2) {
			this->variablePositionPopulation[i] = i;
			this->variablePositionPopulation[j] = j;
		}
		if(i < this->chromosomeSize) {
			this->variablePositionPopulation[i] = i;
		}
}
//
//void Control::shuffleVariablePositionPopulation(SynchronizedUnifGenerator__0__1& unifGen, uint16_t length) {
//	uint16_t randPos = 0;
//	uint16_t popSize = this->variablePositionPopulation.size();
//
//	if(length > popSize) {
//		length = popSize;
//	}
//
//	for(uint16_t i = 0; i < length; ++i) {
//		randPos = i + unifGen() * popSize--;
//		std::swap(this->variablePositionPopulation[randPos], this->variablePositionPopulation[i]);
//	}
//}

std::ostream& operator<<(std::ostream &os, const Control &ctrl) {
	os << "Chromosome size: " << ctrl.chromosomeSize << std::endl
		<< "Population size: " << ctrl.populationSize << std::endl
		<< "Number of generations: " << ctrl.numGenerations << std::endl
		<< "Number of elite chromosomes to keep: " << ctrl.elitism << std::endl
		<< "Number of variables set: " << ctrl.minVariables << " to " << ctrl.maxVariables << std::endl
		<< "Maximum number of tries for mating: " << ctrl.maxMatingTries << std::endl
		<< "Mutation probability: " << ctrl.mutationProbability << std::endl
		<< "Number of threads: " << ctrl.numThreads << std::endl
	<< "Verbosity Level: " << ctrl.verbosity << std::endl;

	return os;
}
