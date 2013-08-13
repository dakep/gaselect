#' Get the found subset(s)
#'
#' Get a list of variable indices/names of the found variable subsets
#'
#' This method is used to get the names or indices of the variables used in specified variable subsets
#'
#' @param object The GenAlg object returned by \code{\link{genAlg}}
#' @param names Should the names of the variables
#' @param indices The indices of the subsets or NULL if all subsets should be returned
#' @export
#' @docType methods
#' @rdname subsets-methods
#' @example examples/subsets.R
setGeneric("subsets", function(object, names = TRUE, indices = NULL) { standardGeneric("subsets"); });

#' @export
#' @rdname subsets-methods
#' @aliases subsets,GenAlg,missing,missing-method
setMethod("subsets", signature(object = "GenAlg", names = "missing", indices = "missing"), function(object, names, indices) {
	subsets(object, TRUE, NULL);
});

#' @export
#' @rdname subsets-methods
#' @aliases subsets,GenAlg,missing,numeric-method
setMethod("subsets", signature(object = "GenAlg", names = "missing", indices = "numeric"), function(object, names, indices) {
	subsets(object, TRUE, indices);
});

#' @export
#' @rdname subsets-methods
#' @aliases subsets,GenAlg,missing,NULL-method
setMethod("subsets", signature(object = "GenAlg", names = "missing", indices = "NULL"), function(object, names, indices) {
	subsets(object, TRUE, NULL);
});

#' @export
#' @rdname subsets-methods
#' @aliases subsets,GenAlg,logical,NULL-method
setMethod("subsets", signature(object = "GenAlg", names = "logical", indices = "NULL"), function(object, names, indices) {
	subsets(object, names, seq_len(ncol(object@subsets)));
});

#' @export
#' @rdname subsets-methods
#' @aliases subsets,GenAlg,logical,missing-method
setMethod("subsets", signature(object = "GenAlg", names = "logical", indices = "missing"), function(object, names, indices) {
	subsets(object, names, NULL);
});

#' @export
#' @rdname subsets-methods
#' @aliases subsets,GenAlg,logical,numeric-method
setMethod("subsets", signature(object = "GenAlg", names = "logical", indices = "numeric"), function(object, names, indices) {
	subs <- object@subsets[ , indices];
	if(names == TRUE && !is.null(colnames(object@covariates))) {
		cnames <- colnames(object@covariates);
		lapply(split(subs, rep(seq_len(length(indices)), each = nrow(subs))), function(s) { cnames[s]; });
	} else {
		lapply(split(subs, rep(seq_len(length(indices)), each = nrow(subs))), function(s) { which(s); });
	}
});
