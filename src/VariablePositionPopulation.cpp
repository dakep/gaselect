//
//  VariablePositionPopulation.cpp
//  GenAlgTest
//
//  Created by David Kepplinger on 12.05.2013.
//
//

#include "VariablePositionPopulation.h"

#include <algorithm>

VariablePositionPopulation::VariablePositionPopulation(const uint16_t size) : size(size), unifGen() {
	uint16_t i = 0, j = 1;
	for(; j < this->size; i += 2, j += 2) {
		this->variablePositionPopulation[i] = i;
		this->variablePositionPopulation[j] = j;
	}
	if(i < this->size) {
		this->variablePositionPopulation[i] = i;
	}
}

VariablePositionPopulation::const_iterator VariablePositionPopulation::shuffle(uint16_t length) {
	uint16_t randPos = 0;
	
	if(length > this->size) {
		throw Rcpp::exception("Can not get a random sample with more elements than the population size", __FILE__, __LINE__);
	}
	
	for(uint16_t i = 0; i < length; ++i) {
		randPos = i + this->unifGen() * (this->size - i);
		std::swap(this->variablePositionPopulation[i], this->variablePositionPopulation[randPos]);
	}
	
	return const_iterator(*this, length);
}

// VariablePositionPopulation::iterator::iterator(const VariablePositionPopulation &obj, const uint16_t length) : obj(obj), length(length) {}

VariablePositionPopulation::const_iterator & VariablePositionPopulation::const_iterator::operator++() {
	if(++this->curPos > this->length) {
		this->curPos = this->obj.variablePositionPopulation.size();
	}
	
	return *this;
}