# genAlg <- function(y, X, control = genAlgControl(popSize = ncol(X) / 10), evalControl = genAlgEvalControl()) {
	# if(chromomsomeSize < 2) {
		# stop("The size of the chromosome must be at least 2");
	# }

	# if(popSize < 2) {
		# stop("The population must at least be 2 chromosomes");
	# }

	# if((popSize %% 2) > 1) {
		# warning("Population size must be a multiple of 2. Using ", popSize - 1, " instead");
		# popSize <- popSize - 1;
	# }

	# if(numGenerations < 2) {
		# stop("The number of generations has to be at least 2");
	# }

	# if(onesRatio <= 0 || onesRatio >= 1) {
		# stop("The ratio of one's must be between 0 and 1");
	# }

	# if(mutationProbability < 0 || mutationProbability > 1) {
		# stop("The mutation probability must be between 0 and 1");
	# }

	# if(!is.function(evalFunction)) {
		# stop("The evaluation function is not callable.");
	# }

	# if(elitism < 0) {
		# warning("Assuming elitism of 0");
	# }

	# if(verbosity > 2) {
		# verbosity <- 2;
	# }

	# invisible(.Call('genAlg', as.integer(chromomsomeSize), as.integer(popSize), as.integer(numGenerations), onesRatio, mutationProbability, as.integer(elitism), as.integer(verbosity), evalFunction, PACKAGE = 'GenAlgPLS'));
# };

simpls <- function(X, Y, ncomp, newX) {
	if(!is.matrix(X) || !is.numeric(X)) {
		stop("X must be a numeric matrix");
	}
	if(!is.matrix(newX) || !is.numeric(newX)) {
		stop("newX must be a numeric matrix");
	}

	if(!(is.matrix(Y) || is.vector(Y)) || !is.numeric(Y)) {
		stop("Y must be a numeric matrix");
	}

	if(is.vector(Y)) {
		Y <- as.matrix(Y);
	}

	ncomp <- as.integer(ncomp);

	if(ncomp < 1 || ncomp >= 2^16) {
		stop("ncomp must be an integer between 1 and ", 2^16);
	}

	invisible(.Call('simpls', X, Y, ncomp, newX, PACKAGE = 'GenAlgPLS'));
};

evalTest <- function(X, Y, numReplications, numSegments) {
	if(!is.matrix(X) || !is.numeric(X)) {
		stop("X must be a numeric matrix");
	}

	if(!((is.matrix(Y) && ncol(Y) == 1) || is.vector(Y)) || !is.numeric(Y)) {
		stop("Y must be a numeric vector or matrix with one column");
	}

	if(is.vector(Y)) {
		Y <- as.matrix(Y);
	}

	numReplications <- as.integer(numReplications);
	numSegments <- as.integer(numSegments);

	if(numReplications < 1 || numReplications > 2^16) {
		stop("numReplications must be between 1 and ", 2^16);
	}
	if(numSegments < 1 || numSegments > 2^16) {
		stop("numSegments must be between 1 and ", 2^16);
	}

	invisible(.Call('evalTest', X, Y, numReplications, numSegments, PACKAGE = 'GenAlgPLS'));
};



genAlg <- function(y, X, control = genAlgControl(populationSize = floor(sqrt(ncol(X))), numGenerations = 100L), evalControl = genAlgEvalControl()) {
	if(!is.numeric(y) || !(is.vector(y) || is.matrix(y) && ncol(y) == 1)) {
		stop("y must be a numeric vector or numeric matrix with 1 column");
	}

	if(!is.numeric(X) || !is.matrix(X)) {
		stop("X must be a numeric matrix");
	}

	if(class(control) != "genAlg-control") {
		stop("control must be an object with class 'genAlg-control'. Use the method 'genAlgControl' to obtain an object of this type");
	}

	if(class(evalControl) != "genAlgEvaluator-control") {
		stop("control must be an object with class 'genAlgEvaluator-control'. Use the method 'genAlgEvalControl' to obtain an object of this type");
	}

	y <- as.matrix(y);

	if(nrow(y) != nrow(X)) {
		stop("Length of y doesn't match number of rows in X");
	}

	ctrlArg <- c(control, evalControl);
	ctrlArg$chromosomeSize = ncol(X);

	if(ctrlArg$minVariables >= ctrlArg$chromosomeSize) {
		stop("The minimum number of variables must be strictly less than the number of available variables");
	}

	if(ctrlArg$maxVariables > ctrlArg$chromosomeSize) {
		stop("The maximum number of variables must be less or equal than the number of available variables");
	}

	if(ctrlArg$useUserSuppliedFunction == TRUE) {
		return(.Call("genAlg", ctrlArg, NULL, NULL, PACKAGE = "GenAlgPLS"));
	} else {
		return(.Call("genAlg", ctrlArg, X, y, PACKAGE = "GenAlgPLS"));
	}
}
