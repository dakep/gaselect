\name{evaluatorPLS}
\alias{evaluatorPLS}
\title{PLS Evaluator}
\usage{
evaluatorPLS(numReplications = 2L, numSegments = 4L, numThreads = NULL,
  method = c("simpls"))
}
\arguments{
  \item{numReplications}{The number of replications used to
  evaluate a variable subset (must be between 1 and 2^16)}

  \item{numSegments}{The number of CV segments used in one
  replication (must be between 1 and 2^16)}

  \item{numThreads}{The maximum number of threads the
  algorithm is allowed to spawn (a value less than 1 or
  NULL means no threads)}

  \item{method}{The PLS method used to fit the PLS model
  (currently only SIMPLS is implemented)}
}
\value{
Returns an S4 object of type
\code{\link{GenAlgPLSEvaluator}}
}
\description{
Create a PLS evaluator for the genetic algorithm
}
\details{
This evaluator class uses PLS with cross-validation (using
\code{numSegments} random segments) to assess the
prediction power of the variable subset. In each of the
\code{numReplications} replications the standard error of
prediction (SEP) is used to quantify the fitness of the
subset. The final fitness is the mean SEP. The larger the
number of replications, the better the estimation of the
SEP but the slower the algorithm (the evaluation step is
done \code{numGenerations} * \code{populationSize} times -
see \code{\link{genAlgControl}}).
}
\examples{
ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorPLS(numReplications = 3L, numSegments = 4L, numThreads = 2L)

\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)
}
}
\seealso{
Other GenAlg.Evaluators: \code{\link{evaluatorLM}};
\code{\link{evaluatorUserFunction}}
}

