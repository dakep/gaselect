#' Evaluator Base Class
#'
#' Virtual base class of all available evaluators
setClass("GenAlgEvaluator", representation(), contains = "VIRTUAL");

#' PLS Evaluator
#'
#' @section Slots:
#' 	\describe{
#' 		\item{\code{numReplications}:}{The number of replications used to evaluate a variable subset.}
#' 		\item{\code{numSegments}:}{The number of CV segments used in one replication.}
#' 		\item{\code{numThreads}:}{The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads).}
#' 		\item{\code{method}:}{The PLS method used to fit the PLS model (currently only SIMPLS is implemented).}
#' 		\item{\code{methodId}:}{The ID of the PLS method used to fit the PLS model.}
#' 	}
#' @rdname GenAlgPLSEvaluator-class
setClass("GenAlgPLSEvaluator", representation(
	numReplications = "integer",
	numSegments = "integer",
	numThreads = "integer",
	method = "character",
	methodId = "integer"
), validity = function(object) {
	errors <- character(0);
	MAXUINT16 <- 2^16; # unsigned 16bit integers are used (uint16_t) in the C++ code

	if(object@numThreads < 0L || object@numThreads > MAXUINT16) {
		errors <- c(errors, paste("The maximum number of threads must be greater than or equal 0 and less than", MAXUINT16));
	}

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
},contains = "GenAlgEvaluator");

#' User Function Evaluator
#'
#' @section Slots:
#' 	\describe{
#' 		\item{\code{evalFunction}:}{The function that is called to evaluate the variable subset.}
#' 		\item{\code{sepFunction}:}{The function that calculates the standard error of prediction for the found subsets.}
#' 	}
#' @rdname GenAlgUserEvaluator-class
setClass("GenAlgUserEvaluator", representation(
	evalFunction = "function",
	sepFunction = "function"
), prototype(
	sepFunction = function(genAlg) {
		warning("Evaluator doesn't support SEP calculation -- using raw fitness");
		return(genAlg@rawFitness);
	}
), contains = "GenAlgEvaluator");

#' LM Evaluator
#'
#' @section Slots:
#' 	\describe{
#' 		\item{\code{statistic}:}{The statistic used to evaluate the fitness.}
#' 		\item{\code{maxCor}:}{If a correlation between the covariates in a subset is higher than this value, the principal
#'								components of the covariates are used to fit the linear model.}
#' 	}
#'
#' @rdname GenAlgLMEvaluator-class
setClass("GenAlgLMEvaluator", representation(
	statistic = "character",
	statisticId = "integer",
	numThreads = "integer"
), prototype(covariatesPC = matrix()), contains = "GenAlgEvaluator",
validity = function(object) {
	errors <- character(0);

	MAXUINT16 <- 2^16; # unsigned 16bit integers are used (uint16_t) in the C++ code

	if(object@numThreads < 0L || object@numThreads > MAXUINT16) {
		errors <- c(errors, paste("The maximum number of threads must be greater than or equal 0 and less than", MAXUINT16));
	}

	if(length(errors) == 0) {
		return(TRUE);
	} else {
		return(errors);
	}
});

#' PLS Evaluator
#'
#' Create a PLS evaluator for the genetic algorithm
#'
#' This evaluator class uses PLS with cross-validation (using \code{numSegments} random segments) to
#' assess the prediction power of the variable subset. In each of the \code{numReplications} replications
#' the standard error of prediction (SEP) is used to quantify the fitness of the subset. The final fitness
#' is the mean SEP. The larger the number of replications, the better the estimation of the SEP but the slower
#' the algorithm (the evaluation step is done \code{numGenerations} * \code{populationSize} times - see \code{\link{genAlgControl}}).
#'
#' @param numReplications The number of replications used to evaluate a variable subset (must be between 1 and 2^16)
#' @param numSegments The number of CV segments used in one replication (must be between 1 and 2^16)
#' @param numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads)
#' @param method The PLS method used to fit the PLS model (currently only SIMPLS is implemented)
#' @export
#' @family GenAlg Evaluators
#' @example examples/genAlg.R
evaluatorPLS <- function(numReplications = 2L, numSegments = 4L, numThreads = NULL, method = c("simpls")) {
	method <- match.arg(method);

	methodId <- switch(method,
		simpls = 0L
	);

	if(is.numeric(numReplications)) {
		numReplications <- as.integer(numReplications);
	}

	if(is.numeric(numSegments)) {
		numSegments <- as.integer(numSegments);
	}

	if(missing(numThreads) || is.null(numThreads)) {
		numThreads <- 1L;
	} else if(is.numeric(numThreads)) {
		numThreads <- as.integer(numThreads);
	}

	return(new("GenAlgPLSEvaluator",
		numReplications = numReplications,
		numSegments = numSegments,
		numThreads = numThreads,
		method = method,
		methodId = methodId
	));
};

#' User Defined Evaluator
#'
#' Create an evaluator that uses a user defined function to evaluate the fitness
#'
#' The user specified function must take a the response vector as first and the covariates matrix as second argument.
#' The function must return a number representing the fitness of the variable subset (the higher the value the fitter the subset)
#' Additionally the user can specify a function that takes a \code{\link{GenAlg}} object and returns
#' the standard error of prediction of the found variable subsets.
#'
#' @param FUN Function used to evaluate the fitness
#' @param sepFUN Function to calculate the SEP of the variable subsets
#' @param ... Additional arguments passed to FUN and sepFUN
#' @export
#' @family GenAlg Evaluators
#' @example examples/evaluatorUserFunction.R
evaluatorUserFunction <- function(FUN, sepFUN = NULL, ...) {
	if(!is.function(FUN)) {
		stop("FUN must be of type `function`");
	};

	evalFunction <- function(y, X) {
		FUN(y, X, ...);
	};

	if(!missing(sepFUN) && is.function(sepFUN)) {
		return(new("GenAlgUserEvaluator",
			evalFunction = evalFunction,
			sepFunction = function(object, genAlg) {
				sepFUN(genAlg, ...);
			}
		));
	} else {
		return(new("GenAlgUserEvaluator",
			evalFunction = evalFunction
		));
	}
};

#' LM Evaluator
#'
#' Create an evaluator that uses a linear model to evaluate the fitness.
#'
#' Different statistics to evaluate the fitness of the variable subset can be given
#'
#' @param statistic The statistic used to evaluate the fitness
#' @param maxCor If the correlation-matrix of the covariates has an entry (absolutely) greater than this value
#'			the principal components are used to calculate the fit
#' @export
#' @family GenAlg Evaluators
#' @example examples/evaluatorLM.R
evaluatorLM <- function(statistic = c("BIC", "AIC", "adjusted.r.squared", "r.squared"), numThreads = NULL, maxCor = NULL) {
	statistic <- match.arg(statistic);

	if(!missing(maxCor) && !is.null(maxCor) && (maxCor < 0 || maxCor > 1)) {
		maxCor <- NULL;
	}

	if(!is.null(maxCor)) {
		FUN <- function(y, X) {
			corrs <- cor(X);
			if(any(abs(corrs[upper.tri(corrs)]) > maxCor)) { #PCA
				X <- prcomp(X)$x;
			}
			m <- lm(y ~ X);
			mss <- sum((m$fitted.values - mean(m$fitted.values)) ^ 2);
			rss <- sum(m$residuals ^ 2);

			switch(statistic,
				adjusted.r.squared = 1 - (1 - (mss / (mss + rss))) * ((nrow(m$qr$qr) - 1) / m$df.residual),
				r.squared = mss / (mss + rss),
				BIC = -BIC(m),
				AIC = -AIC(m)
			)
		};

		sepFUN <- function(genAlg) {
			apply(genAlg@subsets, 2, function(subset) {
				m <- lm(genAlg@response ~ genAlg@covariates[, subset]);
				return(sd(m$residuals));
			});
		};

		return(new("GenAlgUserEvaluator",
			evalFunction = FUN,
			sepFunction = sepFUN
		));

	} else {
		statId = switch(statistic,
			BIC = 0L,
			AIC = 1L,
			adjusted.r.squared = 2L,
			r.squared = 3L
		);

		if(missing(numThreads) || is.null(numThreads)) {
			numThreads <- 1L;
		} else if(is.numeric(numThreads)) {
			numThreads <- as.integer(numThreads);
		}

		return(new("GenAlgLMEvaluator",
			statistic = statistic,
			statisticId = statId,
			numThreads = numThreads
		));
	}
};
