//
//  SingleThreadPopulation.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 16.07.2013.
//
//

#ifndef GenAlgPLS_SingleThreadPopulation_h
#define GenAlgPLS_SingleThreadPopulation_h

#include "config.h"

#include <iostream>
#include <vector>
#include <set>

#include "Chromosome.h"
#include "Evaluator.h"
#include "Control.h"
#include "Population.h"

class SingleThreadPopulation : public Population {
public:
	SingleThreadPopulation(const Control &ctrl, ::Evaluator &evaluator);
	~SingleThreadPopulation() {};

	void run();
private:
};


#endif
