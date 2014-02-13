#' Get the found subset(s)
#'
#' Get a list of variable indices/names of the found variable subsets
#'
#' This method is used to get the names or indices of the variables used in specified variable subsets
#'
#' @param object The GenAlg object returned by \code{\link{genAlg}}
#' @param indices The indices of the subsets or NULL if all subsets should be returned
#' @param names Should the names of the variables
#' @return A logical matrix where each column represents a variable subset
#' @export
#' @include genAlg.R
#' @docType methods
#' @rdname subsets-methods
#' @example examples/subsets.R
setGeneric("subsets", function(object, indices = NULL, names = TRUE) { standardGeneric("subsets"); });

#' @rdname subsets-methods
setMethod("subsets", signature(object = "GenAlg", indices = "missing", names = "missing"), function(object, indices, names) {
	subsets(object, NULL, TRUE);
});

#' @rdname subsets-methods
setMethod("subsets", signature(object = "GenAlg", indices = "numeric", names = "missing"), function(object, indices, names) {
	subsets(object, indices, TRUE);
});

#' @rdname subsets-methods
setMethod("subsets", signature(object = "GenAlg", indices = "NULL", names = "missing"), function(object, indices, names) {
	subsets(object, NULL, TRUE);
});

#' @rdname subsets-methods
setMethod("subsets", signature(object = "GenAlg", indices = "NULL", names = "logical"), function(object, indices, names) {
	subsets(object, seq_len(ncol(object@subsets)), names);
});

#' @rdname subsets-methods
setMethod("subsets", signature(object = "GenAlg", indices = "missing", names = "logical"), function(object, indices, names) {
	subsets(object, NULL, names);
});

#' @rdname subsets-methods
setMethod("subsets", signature(object = "GenAlg", indices = "numeric", names = "logical"), function(object, indices, names) {
	subs <- object@subsets[ , indices, drop = FALSE];
	if(names == TRUE && !is.null(colnames(object@covariates))) {
		cnames <- colnames(object@covariates);
		lapply(split(subs, rep(seq_len(length(indices)), each = nrow(subs))), function(s) { cnames[s]; });
	} else {
		lapply(split(subs, rep(seq_len(length(indices)), each = nrow(subs))), function(s) { which(s); });
	}
});
