#' Evaluate the fitness of variables
#'
#' Evaluate the given data with the given evaluator
#'
#' @param object The GenAlgEvaluator object that is used to evaluate the variables
#' @param X The data matrix used to for fitting the model
#' @param y The response vector
#' @param seed The value to seed the random number generator before evaluating
#' @include Evaluator.R
#' @docType methods
#' @rdname evaluate-methods
setGeneric("evaluate", function(object, X, y, seed) { standardGeneric("evaluate"); });

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "ANY", X = "matrix", y = "numeric", seed = "missing"), function(object, X, y, seed) {
	evaluate(object, X, y, as.integer(sample.int(2^30, 1)));
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "ANY", X = "matrix", y = "numeric", seed = "NULL"), function(object, X, y, seed) {
	evaluate(object, X, y, as.integer(sample.int(2^30, 1)));
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "ANY", X = "matrix", y = "numeric", seed = "numeric"), function(object, X, y, seed) {
	evaluate(object, X, y, as.integer(seed));
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", seed = "integer"), function(object, X, y, seed) {
	ctrlArg <- toCControlList(object);
	ctrlArg$userEvalFunction <- getEvalFun(object, cbind(y, X));
	res <- .Call("evaluate", ctrlArg, as.matrix(X), as.matrix(y), seed, PACKAGE = "GenAlgPLS");
	return(res);
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgUserEvaluator", X = "matrix", y = "numeric", seed = "integer"), function(object, X, y, seed) {
	ctrlArg <- toCControlList(object);
	ctrlArg$userEvalFunction <- getEvalFun(object, cbind(y, X));
	res <- ctrlArg$userEvalFunction(rep(TRUE, ncol(X)));
	return(res);
});
