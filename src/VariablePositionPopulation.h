//
//  VariablePositionPopulation.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 12.05.2013.
//
//

#ifndef GenAlgPLS_VariablePositionPopulation_h
#define GenAlgPLS_VariablePositionPopulation_h

#include "config.h"

#include <iterator>
#include <vector>
#include <RcppArmadillo.h>

#include "RNG.h"

class VariablePositionPopulation;

class VariablePositionPopulation
{
public:
	VariablePositionPopulation(const uint16_t size);

	class const_iterator : public std::iterator<std::input_iterator_tag, uint16_t> {
	public:
		const_iterator(const VariablePositionPopulation &obj, const uint16_t length, const uint16_t shift, const uint16_t pos = 0) : obj(obj), length(length), shift(shift), curPos(pos) {};
		const uint16_t operator*() const { return this->obj.variablePositionPopulation[curPos] + this->shift; }
		const_iterator & operator++();
		bool operator==(const const_iterator & iter) const { return (this->curPos == iter.curPos); }
		bool operator!=(const const_iterator & iter) const { return !(*this == iter); }

	private:
		const VariablePositionPopulation &obj;
		const uint16_t length;
		const uint16_t shift;
		uint16_t curPos;
	};

	typedef VariablePositionPopulation::const_iterator const_iterator;

	/**
	 * Convenience method to shuffle the variable-position population
	 * optionally only the first `length` positions in the population
	 * are random (though all positions are considered)
	 *
	 * The adjustment is added to every index. This is useful if the actual
	 * positions are shifted.
	 */
	const_iterator shuffle(const uint16_t length, const uint16_t shift, RNG& rng);

	const_iterator end() { return const_iterator(*this, 0, 0, this->variablePositionPopulation.size()); };

private:
	const uint16_t size;
	std::vector<uint16_t> variablePositionPopulation;
};


#endif
