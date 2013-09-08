#' Get the evaluation function from a GenAlgUserEvaluator
#'
#' This method returns the correct evaluation function from a GenAlgUserEvaluator
#' that can be used by the C++-code as callback or NULL for any other evaluator
#'
#' @param object The evaluator
#' @param genAlg The \code{\link{GenAlg}} object
#' @docType methods
#' @rdname GenAlgEvaluator-getEvalFun-methods
setGeneric("getEvalFun", function(object, genAlg) { standardGeneric("getEvalFun"); });

#' @rdname GenAlgEvaluator-getEvalFun-methods
#' @aliases getEvalFun,GenAlgUserEvaluator,GenAlg-method
setMethod("getEvalFun", signature(object = "GenAlgUserEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	return(function(varSubset) {
		return(object@evalFunction(genAlg@response, genAlg@covariates[ , varSubset, drop = FALSE]));
	});
});

#' @rdname GenAlgEvaluator-getEvalFun-methods
#' @aliases getEvalFun,GenAlgUserEvaluator,matrix-method
setMethod("getEvalFun", signature(object = "GenAlgUserEvaluator", genAlg = "matrix"), function(object, genAlg) {
	X <- genAlg[ , -1];
	y <- genAlg[ , 1];
	return(function(varSubset) {
		return(object@evalFunction(y, X[ , varSubset, drop = FALSE]));
	});
});

# # @rdname GenAlgEvaluator-getEvalFun-methods
# # @aliases getEvalFun,GenAlgUserEvaluator,GenAlg-method
# setMethod("getEvalFun", signature(object = "GenAlgLMEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
# 	return(function(varSubset) {
# 		return(object@evalFunction(genAlg@response, X));
# 	});
# });

#' @rdname GenAlgEvaluator-getEvalFun-methods
#' @aliases getEvalFun,GenAlgEvaluator,GenAlg-method
setMethod("getEvalFun", signature(object = "GenAlgEvaluator", genAlg = "GenAlg"), function(object, genAlg) {
	return(NULL);
});

#' @rdname GenAlgEvaluator-getEvalFun-methods
#' @aliases getEvalFun,GenAlgEvaluator,matrix-method
setMethod("getEvalFun", signature(object = "GenAlgEvaluator", genAlg = "matrix"), function(object, genAlg) {
	return(NULL);
});
