//
//  ShuffledSet.h
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

#include "RNG.h"

class ShuffledSet
{
public:
	ShuffledSet();
	ShuffledSet(uint16_t size);
//	~ShuffledSet();
	
	class iterator : std::iterator<std::input_iterator_tag, uint16_t> {
	public:
		/**
		 * Attention: The first element in obj.set must already
		 * be a shuffled element!
		 */
		iterator(ShuffledSet &obj, RNG &rng) : obj(obj), rng(rng), pos(0) {
#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
			this->shifted = false;
#endif
		};
		
		const uint16_t operator*() const;
		iterator& operator++();
		bool operator==(const iterator &it) const;
		bool operator!=(const iterator &it) const;

		/**
		 * This method invalidates the shuffle-process
		 * Use only for comparisons (all checks for misuse are disabled!)
		 */
		iterator operator+(const uint16_t &shift);
	private:
		ShuffledSet &obj;
		RNG &rng;
		uint16_t pos;

#ifdef SHUFFLED_SET_CHECK_ITERATOR_STATE
		bool shifted;
#endif
	};

	/**
	 *
	 */
	iterator shuffle(RNG &rng);

	/**
	 * First reset the size of the set to `size`
	 *
	 * @param uint16_t size ... The new size of the set
	 * @param RNG rng ... The random number generator instance
	 * @param bool onlyOne ... If this is TRUE, only one element is guaranteed to be shuffled
	 */
	iterator shuffle(uint16_t size, RNG &rng, bool onlyOne = false);
	
	/**
	 * Reset the set to a sorted stated with `size` elements
	 * (i.e. have values 0, 1, 2, 3, ..., size - 1)
	 *
	 * @param uint16_t size ... The new size of the set
	 */
	void reset(uint16_t size);

	/**
	 * Reset the set to a sorted state with the same number of elements
	 * as before
	 */
	void reset();

private:
	std::vector<uint16_t> set;
};

#endif
