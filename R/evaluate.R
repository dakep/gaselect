#' Evaluate the fitness of variables
#'
#' Evaluate the given data with the given evaluator
#'
#' @param object The GenAlgEvaluator object that is used to evaluate the variables
#' @docType methods
#' @rdname evaluate-methods
setGeneric("evaluate", function(object, X, y, seed = NULL) { standardGeneric("evaluate"); });

#' @rdname evaluate-methods
#' @aliases evaluate,ANY,matrix,numeric,NULL-method
setMethod("evaluate", signature(object = "ANY", X = "matrix", y = "numeric", seed = "NULL"), function(object, X, y, seed) {
	evaluate(object, X, y, as.integer(sample.int(2^30, 1)));
});

#' @rdname evaluate-methods
#' @aliases evaluate,ANY,matrix,numeric,numeric-method
setMethod("evaluate", signature(object = "ANY", X = "matrix", y = "numeric", seed = "numeric"), function(object, X, y, seed) {
	evaluate(object, X, y, as.integer(seed));
});

#' @rdname evaluate-methods
#' @aliases evaluate,GenAlgEvaluator,matrix,numeric,integer-method
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", seed = "integer"), function(object, X, y, seed) {
	ctrlArg <- toCControlList(object);
	ctrlArg$userEvalFunction <- getEvalFun(object, cbind(y, X));
	res <- .Call("evaluate", ctrlArg, as.matrix(X), as.matrix(y), seed, PACKAGE = "GenAlgPLS");
	return(res);
});

#' @rdname evaluate-methods
#' @aliases evaluate,GenAlgUserEvaluator,matrix,numeric,integer-method
setMethod("evaluate", signature(object = "GenAlgUserEvaluator", X = "matrix", y = "numeric", seed = "integer"), function(object, X, y, seed) {
	ctrlArg <- toCControlList(object);
	ctrlArg$userEvalFunction <- getEvalFun(object, cbind(y, X));
	res <- ctrlArg$userEvalFunction(rep(TRUE, ncol(X)));
	return(res);
});
