#' Get the evolution of the fitness
#'
#' Get the fitness of the best / average chromosomes after each generation
#'
#' Returns the progress of the fitness of the best or average chromosome.
#'
#' @param object The \code{\link{GenAlg}} object returned by \code{\link{genAlg}}
#' @param type can be \code{"best"} to return the fitness of the best chromosome for each generation
#' or \code{"avg"} to return the average fitness during each generation
#' @return A vector with the best or average fitness value after each generation
#' @export
#' @include genAlg.R Evaluator.R
#' @docType methods
#' @rdname getFitnessEvolution-methods
setGeneric("getFitnessEvolution", function(object, type = "avg") { standardGeneric("getFitnessEvolution") });

#' @rdname fitness-methods
setMethod("getFitnessEvolution", signature(object = "GenAlg", type = "missing"), function(object, type) {
	return(getFitnessEvolution(object, type));
});
#' @rdname fitness-methods
setMethod("getFitnessEvolution", signature(object = "GenAlg", type = "character"), function(object, type) {
    type <- match.arg(type, c("avg", "best"));
    column <- switch(type, avg = "sum", "best");
    fit <- object@rawFitnessEvolution[ , column, drop = TRUE];

    if (type == "avg") {
        fit <- fit / (object@control@populationSize + object@control@elitism);
    }

	return(trueFitnessVal(object@evaluator, fit));
});


#' Get the fitness of a variable subset
#'
#' Get the internal fitness for all variable subsets
#'
#' This method is used to get the fitness of all variable subsets
#' found by the genetic algorithm.
#'
#' @param object The \code{\link{GenAlg}} object returned by \code{\link{genAlg}}
#' @return A vector with the estimated fitness for each solution
#' @export
#' @include genAlg.R Evaluator.R
#' @docType methods
#' @rdname fitness-methods
#' @example examples/fitness.R
setGeneric("fitness", function(object) { standardGeneric("fitness") });

#' @rdname fitness-methods
setMethod("fitness", signature(object = "GenAlg"), function(object) {
    return(trueFitnessVal(object@evaluator, object@fitness));
});

#' Get the transformed fitness values
#'
#' Transform the given fitness values according tho the GenAlgEvaluator class
#'
#' This method is used to calculate the true fitness given the GenAlgEvaluator class (as they use
#' different internal fitness measures)
#'
#' @param object The used evaluator, an object with type or with a subtype of \code{\link{GenAlgEvaluator}}
#' @param fitness A numeric vector of fitnesses
#' @return A vector with the true fitness values
#' @docType methods
#' @rdname trueFitnessVal-methods
setGeneric("trueFitnessVal", function(object, fitness) { standardGeneric("trueFitnessVal"); });

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgPLSEvaluator", fitness = "numeric"), function(object, fitness) {
	return(switch(object@sepTransformation, log = exp((-fitness)), (-fitness)) / object@numReplications);
});

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgUserEvaluator", fitness = "numeric"), function(object, fitness) {
	return(object@sepFunction(fitness));
});

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgLMEvaluator", fitness = "numeric"), function(object, fitness) {
	return(fitness);
});

#' @rdname trueFitnessVal-methods
setMethod("trueFitnessVal", signature(object = "GenAlgFitEvaluator", fitness = "numeric"), function(object, fitness) {
	return(fitness);
});

