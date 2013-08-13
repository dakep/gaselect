#' Transform the object to a list
#'
#' Get the control list for the C++ procedure genAlgPLS from the object
#'
#' @param object The object
#' @docType methods
#' @rdname toCControlList-methods
setGeneric("toCControlList", function(object) { standardGeneric("toCControlList"); });

#' @rdname toCControlList-methods
#' @aliases toCControlList,GenAlgPLSEvaluator-method
setMethod("toCControlList", signature(object = "GenAlgPLSEvaluator"), function(object) {
	return(list(
		"useUserSuppliedFunction" = FALSE,
		"numReplications" = object@numReplications,
		"numSegments" = object@numSegments,
		"plsMethod" = object@methodId,
		"numThreads" = object@numThreads,
		"userEvalFunction" = function() {NULL;}
	));
});

#' @rdname toCControlList-methods
#' @aliases toCControlList,GenAlgUserEvaluator-method
setMethod("toCControlList", signature(object = "GenAlgUserEvaluator"), function(object) {
	return(list(
		"useUserSuppliedFunction" = TRUE,
		"numReplications" = 0L,
		"numSegments" = 0L,
		"plsMethod" = 0L,
		"numThreads" = 1L,
		"userEvalFunction" = object@evalFunction
	));
});

#' @rdname toCControlList-methods
#' @aliases toCControlList,GenAlgControl-method
setMethod("toCControlList", signature(object = "GenAlgControl"), function(object) {
	return(list(
		"populationSize" = object@populationSize,
		"numGenerations" = object@numGenerations,
		"minVariables" = object@minVariables,
		"maxVariables" = object@maxVariables,
		"maxMatingTries" = object@maxMatingTries,
		"elitism" = object@elitism,
		"mutationProb" = object@mutationProbability,
		"verbosity" = object@verbosity
	));
});
