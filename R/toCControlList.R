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
		"evaluatorClass" = 1,
		"numReplications" = object@numReplications,
		"numSegments" = object@numSegments,
		"plsMethod" = object@methodId,
		"numThreads" = object@numThreads,
		"userEvalFunction" = function() {NULL;},
		"statistic" = 0
	));
});

#' @rdname toCControlList-methods
#' @aliases toCControlList,GenAlgUserEvaluator-method
setMethod("toCControlList", signature(object = "GenAlgUserEvaluator"), function(object) {
	return(list(
		"evaluatorClass" = 0,
		"numReplications" = 0L,
		"numSegments" = 0L,
		"plsMethod" = 0L,
		"numThreads" = 1L,
		"userEvalFunction" = object@evalFunction,
		"statistic" = 0
	));
});

#' @rdname toCControlList-methods
#' @aliases toCControlList,GenAlgLMEvaluator-method
setMethod("toCControlList", signature(object = "GenAlgLMEvaluator"), function(object) {
	return(list(
		"evaluatorClass" = 2,
		"numReplications" = 0L,
		"numSegments" = 0L,
		"plsMethod" = 0L,
		"numThreads" = object@numThreads,
		"userEvalFunction" = function() {NULL;},
		"statistic" = object@statisticId
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
		"crossover" = object@crossoverId,
		"maxDuplicateEliminationTries" = object@maxDuplicateEliminationTries,
		"badSolutionThreshold" = object@badSolutionThreshold,
		"verbosity" = object@verbosity
	));
});
