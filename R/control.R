
genAlgEvalControl <- function(numReplications = 2L, numSegments = 4L, plsMethod = c("simpls"), userEvalFUN = NULL, ...) {
	useUserSuppliedFunction <- (!missing(userEvalFUN) && is.function(userEvalFUN));
	realUserEvalFUN <- NULL;
	if(useUserSuppliedFunction == TRUE) {
		# Proxy additional arguments to user function
		realUserEvalFUN <- function(variables) {
			userEvalFUN(variables, ...);
		} ;
	}

	plsMethod <- match.arg(plsMethod);

	plsMethodNum <- switch(plsMethod,
		simpls = 0L
	);

	ctrl <- list(
		"useUserSuppliedFunction" = useUserSuppliedFunction,
		"numReplications" = as.integer(numReplications),
		"numSegments" = as.integer(numSegments),
		"plsMethod" = plsMethodNum,
		"userEvalFunction" = realUserEvalFUN
	);

	class(ctrl) <- "genAlgEvaluator-control";
	return(ctrl);
};

genAlgControl <- function(populationSize, numGenerations, minVariables, maxVariables, elitism = 10L, mutationProbability = 0.01, verbosity = 0L) {
	MAXUINT16 <- 2^16; # Internal unsigned 16bit integers are used (uint16_t)

	## Type checks:
	if(populationSize < 0L || populationSize > MAXUINT16) {
		stop("The population size must be between 0 and ", MAXUINT16);
	}

	if(numGenerations < 0L || numGenerations > MAXUINT16) {
		stop("The number of generations must be between 0 and ", MAXUINT16);
	}

	if(elitism < 0L || elitism > MAXUINT16) {
		warning("Not storing best solutions from all generations. 'elitism' should be between 0 and ", MAXUINT16);
		elitism <- 0L;
	}

	if(minVariables < 0L || minVariables > MAXUINT16) {
		stop("The minimal number of variables must be between 0 and ", MAXUINT16);
		minVariables <- 0L;
	}

	if(maxVariables < 0L || maxVariables > MAXUINT16) {
		stop("The maximum number of variables must be strictly larger than the minimum number of variables and between 0 and ", MAXUINT16);
		maxVariables <- 0L;
	}

	populationSize  <- as.integer(populationSize);
	numGenerations  <- as.integer(numGenerations);
	elitism <- as.integer(elitism);
	minVariables <- as.integer(minVariables);
	maxVariables <- as.integer(maxVariables);

	## Sanity checks:

	if(minVariables >= maxVariables) {
		stop("The minimal number of variables must be strictly less than the maximum number");
	}

	if(elitism >= populationSize * numGenerations) {
		warning("Requested more elite solutions than possible");
		elitism <- populationSize * numGenerations;
	}

	if(mutationProbability < 0 || mutationProbability >= 1) {
		stop("The mutation probability must be between 0 and 1 (excluded).");
	}

	if(verbosity < 0) {
		verbosity <- 0;
		warning("The verbosity level can not be less than 0 and thus will be set to 0");
	}

	if(verbosity > 3) {
		verbosity <- 3;
		warning("The verbosity level can not be greater than 3 and will thus be set to 3");
	}

	ctrl <- list(
		"populationSize" = populationSize,
		"numGenerations" = numGenerations,
		"minVariables" = minVariables,
		"maxVariables" = maxVariables,
		"elitism" = elitism,
		"mutationProb" = mutationProbability,
		"verbosity" = as.integer(verbosity)
	);

	class(ctrl) <- "genAlg-control";
	return(ctrl);
};