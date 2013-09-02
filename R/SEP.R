#' Get the SEP of a variable subset
#'
#' Calculate the standard error of prediction for all variable subsets
#'
#' This method is used to get the standard error of prediction (SEP) of all variable subsets
#' found by the genetic algorithm.
#' The evaluator may not support the calculation of the SEP. In this case, the raw fitness
#' returned by the evaluator is returned and a warning is issued
#'
#' @param object The GenAlg object returned by \code{\link{genAlg}}
#' @export
#' @docType methods
#' @rdname SEP-methods
#' @example examples/SEP.R
setGeneric("SEP", function(object) { standardGeneric("SEP"); });

#' @export
#' @rdname SEP-methods
#' @aliases SEP,GenAlg-method
setMethod("SEP", signature(object = "GenAlg"), function(object) {
	return(SEPeval(object@evaluator, object));
});

#' Get the SEP of all variable subsets
#'
#' Calculate the standard error of prediction for all variable subsets
#'
#' This method is used to get the standard error of prediction (SEP) of all variable subsets
#' found by the genetic algorithm.
#'
#' @param object The used evaluator
#' @param genAlg The GenAlg object
#' @docType methods
#' @rdname SEPeval-methods
setGeneric("SEPeval", function(object, genAlg) { standardGeneric("SEPeval"); });

#' @rdname SEPeval-methods
#' @aliases SEPeval,GenAlgPLSEvaluator,GenAlg-method
setMethod("SEPeval", signature(object = "GenAlgPLSEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	sumSEP <- (-genAlg@rawFitness);
	return(sumSEP / object@numReplications);
});

#' @rdname SEPeval-methods
#' @aliases SEPeval,GenAlgUserEvaluator,GenAlg-method
setMethod("SEPeval", signature(object = "GenAlgUserEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	return(object@sepFunction(genAlg));
});

setMethod("SEPeval", signature(object = "GenAlgLMEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	return(genAlg@rawFitness);
});