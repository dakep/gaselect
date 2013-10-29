\docType{methods}
\name{SEP}
\alias{SEP}
\alias{SEP,GenAlg-method}
\title{Get the SEP of a variable subset}
\arguments{
  \item{object}{The GenAlg object returned by
  \code{\link{genAlg}}}
}
\description{
  Calculate the standard error of prediction for all
  variable subsets
}
\details{
  This method is used to get the standard error of
  prediction (SEP) of all variable subsets found by the
  genetic algorithm. The evaluator may not support the
  calculation of the SEP. In this case, the raw fitness
  returned by the evaluator is returned and a warning is
  issued
}
\examples{
ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorPLS(numReplications = 3L, numSegments = 4L, numThreads = 2L)

\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)

SEP(result)
}
}

