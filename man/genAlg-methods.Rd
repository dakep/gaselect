\docType{methods}
\name{genAlg}
\alias{genAlg}
\alias{genAlg,numeric,data.frame,GenAlgControl,ANY,ANY-method}
\alias{genAlg,numeric,matrix,GenAlgControl,ANY,ANY-method}
\alias{genAlg,numeric,matrix,GenAlgControl,GenAlgEvaluator,integer-method}
\alias{genAlg,numeric,matrix,GenAlgControl,GenAlgEvaluator,NULL-method}
\alias{genAlg,numeric,matrix,GenAlgControl,GenAlgEvaluator,numeric-method}
\title{Genetic algorithm for variable subset selection}
\arguments{
  \item{y}{The numeric response vector of length n}

  \item{X}{A n x p numeric matrix with all p covariates}

  \item{control}{Options for controlling the genetic
  algorithm. See \code{\link{genAlgControl}} for details.}

  \item{evaluator}{The evaluator used to evaluate the
  fitness of a variable subset. See
  \code{\link{evaluatorPLS}}, \code{\link{evaluatorLM}} or
  \code{\link{evaluatorUserFunction}} for details.}

  \item{seed}{Integer with the seed for the random number
  generator or NULL to automatically seed the RNG}
}
\description{
  A genetic algorithm to find "good" variable subsets based
  on internal PLS evaluation or a user specified evaluation
  function
}
\details{
  The GA generates an initial "population" of
  \code{populationSize} chromosomes where each initial
  chromosome has a random number of randomly selected
  variables. The fitness of every chromosome is evaluated
  by the specified evaluator. The default built-in PLS
  evaluator (see \code{\link{evaluatorPLS}}) is the
  preferred evaluator. Chromosomes with higher fitness have
  higher probability of mating with another chromosome.
  \code{populationSize / 2} couples each create 2 children.
  The children are created by randomly mixing the parents'
  variables. These children make up the new generation and
  are again selected for mating based on their fitness. A
  total of \code{numGenerations} generations are built this
  way. The algorithm returns the last generation as well as
  the best \code{elitism} chromosomes from all generations.
}
\examples{
ctrl <- genAlgControl(populationSize = 200, numGenerations = 30,
	minVariables = 5, maxVariables = 15)
evaluator <- evaluatorPLS(numReplications = 3L, numSegments = 4L, numThreads = 2L)

\dontrun{
result <- genAlg(y, x, control = ctrl, evaluator = evaluator)
}
}

