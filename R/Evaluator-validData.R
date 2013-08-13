#' Check if the data is valid for the evaluator
#'
#' This method checks if the covariates matrix is valid for the evaluator
#'
#' @param object The evaluator
#' @param X The covariates matrix
#' @usage validData(evaluator, X)
#' @docType methods
#' @rdname GenAlgEvaluator-validData-methods
#' @examples
#' evaluator <- evaluatorLM()
#' validData(evaluator, matrix(0, ncol = 100, nrow = 50); # not valid because more variables than observations
setGeneric("validData", function(object, X) { standardGeneric("validData"); });

#' @rdname GenAlgEvaluator-validData-methods
#' @aliases validData,GenAlgPLSEvaluator,ANY-method
setMethod("validData", signature(object = "GenAlgPLSEvaluator", X = "ANY"), function(object, X) {
	if(is.matrix(X) && is.numeric(X)) {
		return(TRUE);
	} else {
		return("The covariates have to be numerical");
	}
});

#' @rdname GenAlgEvaluator-validData-methods
#' @aliases validData,GenAlgLMEvaluator,ANY-method
setMethod("validData", signature(object = "GenAlgLMEvaluator", X = "ANY"), function(object, X) {
	if((is.matrix(X) || is.data.frame(X)) && (ncol(X) < nrow(X))) {
		return(TRUE);
	} else {
		return("It is not possible to use a linear model if there are more observations than covariates.");
	}
});

#' @rdname GenAlgEvaluator-validData-methods
#' @aliases validData,GenAlgEvaluator,ANY-method
setMethod("validData", signature(object = "GenAlgEvaluator", X = "ANY"), function(object, X) {
	return(TRUE);
});
