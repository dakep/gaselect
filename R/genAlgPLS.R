#' Genetic algorithm for variable subset selection
#'
#' A genetic algorithm to find "good" variable subsets based on internal PLS evaluation or a user specified
#' evaluation function
#'
#' The GA generates an initial "population" of \code{populationSize} chromosomes where each initial
#' chromosome has a random number of randomly selected variables. The fitness of every chromosome is evaluated
#' using the built-in PLS method to evaluate the standard error of prediction (SEP) or a user defined evaluation function (see \code{\link{genAlgEvalControl}}).
#' Chromosomes with higher fitness have higher probability of mating with another chromosome. \code{populationSize / 2} couples each create
#' 2 children. The children are created by randomly mixing the parents' variables. These children make up the new generation and are again
#' selected for mating based on their fitness. A total of \code{numGenerations} generations are built this way.
#' The algorithm returns the last generation as well as the best \code{elitism} chromosomes from all generations.
#'
#' @param y The numeric response vector of length n
#' @param X A n x p numeric matrix with all p covariates
#' @param control Options for controlling the genetic algorithm. See \code{\link{genAlgControl}} for details.
#' @param evalControl Options for controlling the evaluation step. See \code{\link{genAlgEvalControl}} for details.
#' @return Returns a list with elements \code{subsets} and \code{fitness}.
#'		\item{\code{subsets}}{logical matrix with one variable subset per column. The columns are ordered according to their fitness (first column contains the fittest variable-subset)}
#' 		\item{\code{fitness}}{numeric vector with the fitness of the corresponding variable subset}
#' @export
#' @useDynLib GenAlgPLS
genAlgPLS <- function(y, X, control = genAlgControl(populationSize = floor(sqrt(ncol(X))), numGenerations = 100L), evalControl = genAlgEvalControl()) {
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
		return(.Call("genAlgPLS", ctrlArg, NULL, NULL, PACKAGE = "GenAlgPLS"));
	} else {
		return(.Call("genAlgPLS", ctrlArg, X, y, PACKAGE = "GenAlgPLS"));
	}
}

# simpls <- function(X, Y, ncomp, newX) {
# 	if(!is.matrix(X) || !is.numeric(X)) {
# 		stop("X must be a numeric matrix");
# 	}
# 	if(!is.matrix(newX) || !is.numeric(newX)) {
# 		stop("newX must be a numeric matrix");
# 	}
#
# 	if(!(is.matrix(Y) || is.vector(Y)) || !is.numeric(Y)) {
# 		stop("Y must be a numeric matrix");
# 	}
#
# 	if(is.vector(Y)) {
# 		Y <- as.matrix(Y);
# 	}
#
# 	ncomp <- as.integer(ncomp);
#
# 	if(ncomp < 1 || ncomp >= 2^16) {
# 		stop("ncomp must be an integer between 1 and ", 2^16);
# 	}
#
# 	invisible(.Call('simpls', X, Y, ncomp, newX, PACKAGE = 'GenAlgPLS'));
# };
#
# evalTest <- function(X, Y, numReplications, numSegments) {
# 	if(!is.matrix(X) || !is.numeric(X)) {
# 		stop("X must be a numeric matrix");
# 	}
#
# 	if(!((is.matrix(Y) && ncol(Y) == 1) || is.vector(Y)) || !is.numeric(Y)) {
# 		stop("Y must be a numeric vector or matrix with one column");
# 	}
#
# 	if(is.vector(Y)) {
# 		Y <- as.matrix(Y);
# 	}
#
# 	numReplications <- as.integer(numReplications);
# 	numSegments <- as.integer(numSegments);
#
# 	if(numReplications < 1 || numReplications > 2^16) {
# 		stop("numReplications must be between 1 and ", 2^16);
# 	}
# 	if(numSegments < 1 || numSegments > 2^16) {
# 		stop("numSegments must be between 1 and ", 2^16);
# 	}
#
# 	invisible(.Call('evalTest', X, Y, numReplications, numSegments, PACKAGE = 'GenAlgPLS'));
# };
