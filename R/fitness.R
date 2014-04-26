#' Get the fitness of a variable subset
#'
#' Get the internal fitness for all variable subsets
#'
#' This method is used to get the fitness of all variable subsets
#' found by the genetic algorithm.
#'
#' @param object The \code{\link{GenAlg}} object returned by \code{\link{genAlg}}
#' @return A vector with the estimated fitness for each solution
#' @export
#' @include genAlg.R Evaluator.R
#' @docType methods
#' @rdname fitness-methods
#' @example examples/fitness.R
setGeneric("fitness", function(object) { standardGeneric("fitness") });

#' @rdname fitness-methods
setMethod("fitness", signature(object = "GenAlg"), function(object) {
	return(fitnessEval(object@evaluator, object));
});

#' Get the fitness of a variable subset
#'
#' Get the internal fitness for all variable subsets
#'
#' This method is used to get the fitness of all variable subsets
#' found by the genetic algorithm.
#'
#' @param object The used evaluator, i.e. an object of type \code{\link{GenAlgEvaluator}}
#' @param genAlg An object of type \code{\link{GenAlg}}
#' @return A vector with the estimated fitness for each solution
#' @docType methods
#' @rdname fitnessEval-methods
setGeneric("fitnessEval", function(object, genAlg) { standardGeneric("fitnessEval"); });

#' @rdname fitnessEval-methods
setMethod("fitnessEval", signature(object = "GenAlgPLSEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	sumSEP <- (-genAlg@rawFitness);
	return(sumSEP / object@numReplications);
});

#' @rdname fitnessEval-methods
setMethod("fitnessEval", signature(object = "GenAlgUserEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	return(object@sepFunction(genAlg));
});

#' @rdname fitnessEval-methods
setMethod("fitnessEval", signature(object = "GenAlgLMEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	return(genAlg@rawFitness);
});
