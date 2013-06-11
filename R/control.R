#' Set control arguments for the evaluation-step in the genetic algorithm
#'
#' Controls the built-in evaluation method or tells the GA to use a user-specified evaluation method.
#'
#' The evaluation of variable subsets is crucial for the performance of the genAlgPLS method.
#' The user can use the built-in evaluation method or an user-supplied function.
#'
#' The built-in method uses PLS with cross-validation (using \code{numSegments} random segments) to
#' assess the prediction power of the variable subset. In each of the \code{numReplications} replications
#' the standard error of prediction (SEP) is used to quantify the fitness of the subset. The final fitness
#' is the mean SEP. The larger the number of replications, the better the estimation of the SEP but the slower
#' the algorithm (the evaluation step is done \code{numGenerations} * \code{populationSize} times - see \code{\link{genAlgControl}}).
#'
#' The user specified function must take a logical vector as its first argument. This logical vector
#' specifies the variables in the subset. The function must return a number representing the fitness
#' of the variable subset (the higher the value the fitter the subset)
#'
#' @param numReplications The number of replications used to evaluate a variable subset (must be between 1 and 2^16)
#' @param numSegments The number of CV segments used in one replication (must be between 1 and 2^16)
#' @param plsMethod The PLS method used to fit the PLS model (currently only SIMPLS implemented)
#' @param userEvalFUN Alternatively the user can specify an evaluation function (see Details)
#' @param ... Additional arguments passed to the user specified evaluation function
#' @export
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

#' Set control arguments for the genetic algorithm
#'
#' The population must be large enough to allow the algorithm to explore the whole solution space. If
#' the initial population is not diverse enough, the chance to find the global optimum is very small.
#' Thus the more variables to choose from, the larger the population has to be.
#'
#' The initial population is generated randomly. Every chromosome uses between \code{minVariables} and
#' \code{maxVariables} (uniformly distributed).
#'
#' If the mutation probability (\code{mutationProbability} is greater than 0, a random number of
#' variables is added/removed according to a truncated geometric distribution to each offspring-chromosome.
#' The resulting distribution of the total number of variables in the subset is not uniform anymore, but almost (the smaller the
#' mutation probability, the more "uniform" the distribution). This should not be a problem for most
#' applications.
#'
#' Elitism is a method of enhancing the GA by keeping track of very good solutions. The parameter \code{elitism}
#' specifies how many "very good" solutions should be kept.
#'
#' @param populationSize The number of "chromosomes" in the population (between 1 and 2^16)
#' @param numGenerations The number of generations to produce (between 1 and 2^16)
#' @param minVariables The minimum number of variables in the variable subset (between 0 and p - 1 where p is the total number of variables)
#' @param maxVariables The maximum number of variables in the variable subset (between 1 and p, and greater than \code{minVariables})
#' @param elitism The number of absolute best chromosomes to keep across all generations (between 1 and min(\code{populationSize} * \code{numGenerations}, 2^16))
#' @param mutationProbability The probability of mutation (between 0 and 1)
#' @param verbosity The level of verbosity. 0 means no output at all, 3 is very verbose.
#' @export
#' @seealso See \code{\link{genAlgEvalControl}} for controlling the evaluation step during the GA
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
