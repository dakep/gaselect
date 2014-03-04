#' Evaluate the fitness of variables
#'
#' Evaluate the given data with the given evaluator
#'
#' @param object The GenAlgEvaluator object that is used to evaluate the variables
#' @param X The data matrix used to for fitting the model
#' @param y The response vector
#' @param subsets The logical matrix where a column stands for one subset to evaluate
#' @param seed The value to seed the random number generator before evaluating
#' @include Evaluator.R
#' @docType methods
#' @rdname evaluate-methods
setGeneric("evaluate", function(object, X, y, subsets, seed) { standardGeneric("evaluate"); });

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "ANY", seed = "missing"), function(object, X, y, subsets, seed) {
	evaluate(object, X, y, subsets, as.integer(sample.int(2^30, 1)));
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "ANY", seed = "NULL"), function(object, X, y, subsets, seed) {
	evaluate(object, X, y, subsets, as.integer(sample.int(2^30, 1)));
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "ANY", seed = "numeric"), function(object, X, y, subsets, seed) {
	evaluate(object, X, y, subsets, as.integer(seed));
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "logical", seed = "integer"), function(object, X, y, subsets, seed) {
	evaluate(object, X, y, as.matrix(subsets), seed);
});

#' @rdname evaluate-methods
setMethod("evaluate", signature(object = "GenAlgEvaluator", X = "matrix", y = "numeric", subsets = "matrix", seed = "integer"), function(object, X, y, subsets, seed) {
	if(!is.logical(subsets)) {
		stop("subsets must be logical.");
	}

	if(!is.numeric(X)) {
		stop("X must be numeric.");
	}

	if(nrow(subsets) != ncol(X)) {
		stop("The number of rows of subsets must match the number of columns of X.");
	}

	ctrlArg <- toCControlList(object);
	ctrlArg$userEvalFunction <- getEvalFun(object, cbind(y, X));
	res <- .Call("evaluate", ctrlArg, as.matrix(X), as.matrix(y), subsets, seed, PACKAGE = "GenAlgPLS");

	if(class(object) == "GenAlgPLSEvaluator") {
		res <- -res / object@numReplications;
	}

	return(res);
});
