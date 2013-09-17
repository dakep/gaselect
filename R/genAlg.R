#' Result of a genetic algorithm run
#'
#' Return object of a run of the genetic algorithm genAlg
#' @section Slots:
#' 	\describe{
#' 		\item{\code{subsets}:}{Logical matrix with one variable subset per column. The columns are ordered according to their fitness (first column contains the fittest variable-subset).}
#' 		\item{\code{rawFitness}:}{Numeric vector with the raw fitness of the corresponding variable subset returned by the evaluator.}
#' 		\item{\code{response}:}{The original response vector.}
#' 		\item{\code{covariates}:}{The original covariates matrix.}
#' 		\item{\code{evaluator}:}{The evaluator used in the genetic algorithm.}
#' 		\item{\code{control}:}{The control object.}
#' 	}
#' @rdname GenAlg-class
setClass("GenAlg", representation(
	subsets = "matrix",
	rawFitness = "numeric",
	response = "numeric",
	covariates = "matrix",
	evaluator = "GenAlgEvaluator",
	control = "GenAlgControl"
), prototype(
	subsets = matrix(),
	rawFitness = NA_real_
), validity = function(object) {
	errors <- character(0);
	if(!is.numeric(object@response) || !(is.vector(object@response) || is.matrix(object@response) && ncol(object@response) == 1)) {
		errors <- c(errors, "The response variable must be a numeric vector");
	}

	if(!is.numeric(object@covariates) || !is.matrix(object@covariates)) {
		errors <- c(errors, "The covariates must be a numerical matrix");
	}

	if(length(object@response) != nrow(object@covariates)) {
		errors <- c(errors, "The response and the covariates must have the same number of observations");
	}

	if(object@control@minVariables >= ncol(object@covariates)) {
		errors <- c(errors, "The minimum number of variables must be strictly less than the number of available variables");
	}

	if(object@control@maxVariables > ncol(object@covariates)) {
		errors <- c(errors, "The maximum number of variables must be less or equal than the number of available variables");
	}

	dataErrors <- validData(object@evaluator, object)
	if(!is.logical(dataErrors)) {
		errors <- c(errors, dataErrors);
	} else if(dataErrors == FALSE) {
		errors <- c(errors, "The evaluator can not handle this kind of data");
	}

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
});

#' Genetic algorithm for variable subset selection
#'
#' A genetic algorithm to find "good" variable subsets based on internal PLS evaluation or a user specified
#' evaluation function
#'
#' The GA generates an initial "population" of \code{populationSize} chromosomes where each initial
#' chromosome has a random number of randomly selected variables. The fitness of every chromosome is evaluated by
#' the specified evaluator. The default built-in PLS evaluator (see \code{\link{evaluatorPLS}}) is the preferred
#' evaluator.
#' Chromosomes with higher fitness have higher probability of mating with another chromosome. \code{populationSize / 2} couples each create
#' 2 children. The children are created by randomly mixing the parents' variables. These children make up the new generation and are again
#' selected for mating based on their fitness. A total of \code{numGenerations} generations are built this way.
#' The algorithm returns the last generation as well as the best \code{elitism} chromosomes from all generations.
#'
#' @param y The numeric response vector of length n
#' @param X A n x p numeric matrix with all p covariates
#' @param control Options for controlling the genetic algorithm. See \code{\link{genAlgControl}} for details.
#' @param evaluator The evaluator used to evaluate the fitness of a variable subset. See \code{\link{evaluatorPLS}}, \code{\link{evaluatorLM}} or \code{\link{evaluatorUserFunction}} for details.
#' @param seed Integer with the seed for the random number generator or NULL to automatically seed the RNG
#' @export
#' @rdname genAlg-methods
#' @docType methods
#' @useDynLib GenAlgPLS
#' @example examples/genAlg.R
setGeneric("genAlg", function(y, X, control, evaluator, seed) { standardGeneric("genAlg") });

#' @rdname genAlg-methods
#' @aliases genAlg,numeric,matrix,GenAlgControl,missing,missing-method
setMethod("genAlg", signature(y = "numeric", X = "matrix", control = "GenAlgControl", evaluator = "missing", seed = "missing"),
function(y, X, control, evaluator) {
	genAlg(y, X, control, evaluatorPLS(), as.integer(sample.int(2^30, 1)));
});

#' @rdname genAlg-methods
#' @aliases genAlg,numeric,matrix,GenAlgControl,missing,ANY-method
setMethod("genAlg", signature(y = "numeric", X = "matrix", control = "GenAlgControl", evaluator = "missing", seed = "ANY"),
function(y, X, control, evaluator, seed) {
	genAlg(y, X, control, evaluatorPLS(), seed);
});

#' @rdname genAlg-methods
#' @aliases genAlg,numeric,matrix,GenAlgControl,GenAlgEvaluator,NULL-method
setMethod("genAlg", signature(y = "numeric", X = "matrix", control = "GenAlgControl", evaluator = "GenAlgEvaluator", seed = "NULL"),
function(y, X, control, evaluator = evaluatorPLS(), seed) {
	genAlg(y, X, control, evaluator, as.integer(sample.int(2^30, 1)));
});
#' @rdname genAlg-methods
#' @aliases genAlg,numeric,matrix,GenAlgControl,GenAlgEvaluator,NULL-method
setMethod("genAlg", signature(y = "numeric", X = "matrix", control = "GenAlgControl", evaluator = "GenAlgEvaluator", seed = "missing"),
function(y, X, control, evaluator = evaluatorPLS(), seed) {
	genAlg(y, X, control, evaluator, as.integer(sample.int(2^30, 1)));
});

#' @rdname genAlg-methods
#' @aliases genAlg,numeric,matrix,GenAlgControl,GenAlgEvaluator,numeric-method
setMethod("genAlg", signature(y = "numeric", X = "matrix", control = "GenAlgControl", evaluator = "GenAlgEvaluator", seed = "numeric"),
function(y, X, control, evaluator = evaluatorPLS(), seed) {
	genAlg(y, X, control, evaluator, as.integer(seed));
});

#' @rdname genAlg-methods
#' @aliases genAlg,numeric,matrix,GenAlgControl,GenAlgEvaluator,integer-method
setMethod("genAlg", signature(y = "numeric", X = "matrix", control = "GenAlgControl", evaluator = "GenAlgEvaluator", seed = "integer"),
function(y, X, control, evaluator = evaluatorPLS(), seed) {
	ret <- new("GenAlg",
		response = y,
		covariates = X,
		evaluator = evaluator,
		control = control
	);

	ctrlArg <- c(toCControlList(ret@control), toCControlList(ret@evaluator));
	ctrlArg$chromosomeSize = ncol(ret@covariates);

	ctrlArg$userEvalFunction <- getEvalFun(ret@evaluator, ret);

	if(ctrlArg$evaluatorClass == 0) {
		res <- .Call("genAlgPLS", ctrlArg, NULL, NULL, seed, PACKAGE = "GenAlgPLS");
	} else {
		res <- .Call("genAlgPLS", ctrlArg, ret@covariates, as.matrix(ret@response), seed, PACKAGE = "GenAlgPLS");
	}

	ret@subsets <- res$subsets;
	ret@rawFitness <- res$fitness;

	return(ret);
});
