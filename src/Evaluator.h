//
//  Evaluator.h
//  GenAlgPLS
//
//  Created by David Kepplinger on 15.04.2013.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//

#ifndef GenAlgPLS_Evaluator_h
#define GenAlgPLS_Evaluator_h

#include "config.h"
#include "Control.h"
#include "Chromosome.h"

class Evaluator {
public:
	Evaluator(const VerbosityLevel verbosity) : verbosity(verbosity) {}
	virtual ~Evaluator() {};

	virtual double evaluate(Chromosome &ch) const = 0;
	
	virtual Evaluator* clone() const = 0;
protected:
	const VerbosityLevel verbosity;
};


#endif
