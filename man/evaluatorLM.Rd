\name{evaluatorLM}
\alias{evaluatorLM}
\title{LM Evaluator}
\usage{
  evaluatorLM(statistic = c("adjusted.r.squared", "r.squared", "BIC", "AIC"),
    maxCor = 0.95)
}
\arguments{
  \item{statistic}{The statistic used to evaluate the
  fitness}

  \item{maxCor}{If the correlation-matrix of the covariates
  has an entry (absolutely) greater than this value the
  principal components are used to calculate the fit}
}
\description{
  Create an evaluator that uses a linear model to evaluate
  the fitness.
}
\details{
  Different statistics to evaluate the fitness of the
  variable subset can be given
}
\examples{
ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorLM(statistic = "adjusted.r.squared", maxCor = 0.8)

\\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)
}
}
\seealso{
  Other GenAlg.Evaluators: \code{\link{evaluatorPLS}},
  \code{\link{evaluatorUserFunction}}
}

