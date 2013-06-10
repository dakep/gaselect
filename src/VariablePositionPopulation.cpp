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
	this->variablePositionPopulation.reserve(this->size);
	for(; j < this->size; i += 2, j += 2) {
		this->variablePositionPopulation.push_back(i);
		this->variablePositionPopulation.push_back(j);
	}
	if(i < this->size) {
		this->variablePositionPopulation.push_back(i);
	}
}

VariablePositionPopulation::const_iterator VariablePositionPopulation::shuffle(const uint16_t length, const uint16_t shift) {
	uint16_t randPos = 0;

	if(length > this->size) {
		throw Rcpp::exception("Can not get a random sample with more elements than elements in the population", __FILE__, __LINE__);
	}

	for(uint16_t i = 0; i < length; ++i) {
		randPos = i + this->unifGen() * (this->size - i);
		std::swap(this->variablePositionPopulation[i], this->variablePositionPopulation[randPos]);
	}

	return const_iterator(*this, length, shift);
}

VariablePositionPopulation::const_iterator & VariablePositionPopulation::const_iterator::operator++() {
	if(++this->curPos >= this->length) {
		this->curPos = this->obj.variablePositionPopulation.size();
	}

	return *this;
}