\name{evaluatorUserFunction}
\alias{evaluatorUserFunction}
\title{User Defined Evaluator}
\usage{
  evaluatorUserFunction(FUN, sepFUN = NULL, ...)
}
\arguments{
  \item{FUN}{Function used to evaluate the fitness}

  \item{sepFUN}{Function to calculate the SEP of the
  variable subsets}

  \item{...}{Additional arguments passed to FUN and sepFUN}
}
\description{
  Create an evaluator that uses a user defined function to
  evaluate the fitness
}
\details{
  The user specified function must take a the response
  vector as first and the covariates matrix as second
  argument. The function must return a number representing
  the fitness of the variable subset (the higher the value
  the fitter the subset) Additionally the user can specify
  a function that takes a \code{\link{GenAlg}} object and
  returns the standard error of prediction of the found
  variable subsets.
}
\examples{
ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)

# Use the BIC of a linear model to evaluate the fitness of a variable subset
evalFun <- function(y, X) {
		return(BIC(lm(y ~ X)));
}

# Dummy function that returns the residuals standard deviation and not the SEP
sepFUN <- function(genAlg) {
		return(apply(genAlg@subsets, 2, function(subset) {
			m <- lm(genAlg@response ~ genAlg@covariates[, subset]);
			return(sd(m$residuals));
		}));
}
evaluator <- evaluatorUserFunction(FUN = evalFun, sepFUN = sepFUN)
\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)
}
}
\seealso{
  Other GenAlg.Evaluators: \code{\link{evaluatorLM}},
  \code{\link{evaluatorPLS}}
}

