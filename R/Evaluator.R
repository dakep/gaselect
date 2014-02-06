#' Evaluator Base Class
#'
#' Virtual base class of all available evaluators
#' @aliases GenAlgEvaluator
#' @rdname GenAlgEvaluator-class
setClass("GenAlgEvaluator", representation(), contains = "VIRTUAL");

#' PLS Evaluator
#' @slot numReplications The number of replications used to evaluate a variable subset.
#' @slot numSegments The number of CV segments used in one replication.
#' @slot numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads).
#' @slot method The PLS method used to fit the PLS model (currently only SIMPLS is implemented).
#' @slot methodId The ID of the PLS method used to fit the PLS model.
#' @aliases GenAlgPLSEvaluator
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
#' @slot evalFunction The function that is called to evaluate the variable subset.
#' @slot sepFunction The function that calculates the standard error of prediction for the found subsets.
#' @aliases GenAlgUserEvaluator
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
#' @slot statistic The statistic used to evaluate the fitness.
#' @slot statisticId The (internal) numeric ID of the statistic
#' @slot numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads).
#' @aliases GenAlgLMEvaluator
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
#' @return Returns an S4 object of type \code{\link{GenAlgPLSEvaluator}}
#' @export
#' @family GenAlg Evaluators
#' @example examples/genAlg.R
#' @rdname GenAlgPLSEvaluator-constructor
evaluatorPLS <- function(numReplications = 30L, numSegments = 5L, numThreads = NULL, method = c("simpls")) {
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
#' @return Returns an S4 object of type \code{\link{GenAlgUserEvaluator}}
#' @export
#' @family GenAlg Evaluators
#' @example examples/evaluatorUserFunction.R
#' @rdname GenAlgUserEvaluator-constructor
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
#' Different statistics to evaluate the fitness of the variable subset can be given. If a maximum
#' absolute correlation is given the algorithm will be very slow (as the C++ implementation can not
#' be used anymore) and multithreading is not available.
#'
#' @param statistic The statistic used to evaluate the fitness
#' @param numThreads The maximum number of threads the algorithm is allowed to spawn (a value less than 1 or NULL means no threads)
#' @param maxCor If the correlation-matrix of the covariates has an entry (absolutely) greater than this value
#'			the principal components are used to calculate the fit
#' @return Returns an S4 object of type \code{\link{GenAlgLMEvaluator}}
#' @export
#' @family GenAlg Evaluators
#' @example examples/evaluatorLM.R
#' @rdname GenAlgLMEvaluator-constructor
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
