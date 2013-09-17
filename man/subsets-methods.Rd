\docType{methods}
\name{subsets}
\alias{subsets}
\alias{subsets,GenAlg,NULL,logical-method}
\alias{subsets,GenAlg,NULL,missing-method}
\alias{subsets,GenAlg,missing,logical-method}
\alias{subsets,GenAlg,missing,missing-method}
\alias{subsets,GenAlg,numeric,logical-method}
\alias{subsets,GenAlg,numeric,missing-method}
\title{Get the found subset(s)}
\arguments{
  \item{object}{The GenAlg object returned by
  \code{\link{genAlg}}}

  \item{indices}{The indices of the subsets or NULL if all
  subsets should be returned}

  \item{names}{Should the names of the variables}
}
\description{
  Get a list of variable indices/names of the found
  variable subsets
}
\details{
  This method is used to get the names or indices of the
  variables used in specified variable subsets
}
\examples{
ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorPLS(numReplications = 3L, numSegments = 4L, numThreads = 2L)

\\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)
subsets(result, names = TRUE, indices = 1:5) # best 5 variable subsets as a list of names
result@subsets[ , 1:5] # best 5 variable subsets as a logical matrix with the subsets in the columns
}
}

