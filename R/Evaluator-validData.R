#' Check if the data is valid for the evaluator
#'
#' This method checks if the covariates matrix is valid for the evaluator
#'
#' @param object The evaluator
#' @param genAlg The GenAlg object the evaluator is used in
#' @docType methods
#' @rdname GenAlgEvaluator-validData-methods
setGeneric("validData", function(object, genAlg) { standardGeneric("validData"); });

#' @rdname GenAlgEvaluator-validData-methods
#' @aliases validData,GenAlgPLSEvaluator,GenAlg-method
setMethod("validData", signature(object = "GenAlgPLSEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	if(is.numeric(genAlg@covariates)) {
		return(TRUE);
	} else {
		return("The covariates have to be numerical");
	}
});

#' @rdname GenAlgEvaluator-validData-methods
#' @aliases validData,GenAlgLMEvaluator,GenAlg-method
setMethod("validData", signature(object = "GenAlgLMEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	if(genAlg@control@maxVariables < nrow(genAlg@covariates)) {
		return(TRUE);
	} else {
		return("It is not possible to use a linear model if maxVariables is greater than the number of observations.");
	}
});

#' @rdname GenAlgEvaluator-validData-methods
#' @aliases validData,GenAlgEvaluator,GenAlg-method
setMethod("validData", signature(object = "GenAlgEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	return(TRUE);
});
