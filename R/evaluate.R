#' Evaluate the fitness of variables
#'
#' Evaluate the given data with the given evaluator
#'
#' @param object The GenAlgEvaluator object that is used to evaluate the variables
#' @docType methods
#' @rdname evaluate-methods
setGeneric("evaluate", function(object, X, y) { standardGeneric("evaluate"); });

#' @rdname evaluate-methods
#' @aliases evaluate,GenAlgEvaluator,matrix,numeric-method
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric"), function(object, X, y) {
	ctrlArg <- toCControlList(object);
	ctrlArg$userEvalFunction <- getEvalFun(object, cbind(y, X));
	res <- .Call("evaluate", ctrlArg, as.matrix(X), as.matrix(y), PACKAGE = "GenAlgPLS");
	return(res);
});
