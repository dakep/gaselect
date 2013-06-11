\name{genAlgEvalControl}
\alias{genAlgEvalControl}
\title{Set control arguments for the evaluation-step in the genetic algorithm}
\usage{
  genAlgEvalControl(numReplications = 2L, numSegments = 4L,
    plsMethod = c("simpls"), userEvalFUN = NULL, ...)
}
\arguments{
  \item{numReplications}{The number of replications used to
  evaluate a variable subset (must be between 1 and 2^16)}

  \item{numSegments}{The number of CV segments used in one
  replication (must be between 1 and 2^16)}

  \item{plsMethod}{The PLS method used to fit the PLS model
  (currently only SIMPLS implemented)}

  \item{userEvalFUN}{Alternatively the user can specify an
  evaluation function (see Details)}

  \item{...}{Additional arguments passed to the user
  specified evaluation function}
}
\description{
  Controls the built-in evaluation method or tells the GA
  to use a user-specified evaluation method.
}
\details{
  The evaluation of variable subsets is crucial for the
  performance of the genAlgPLS method. The user can use the
  built-in evaluation method or an user-supplied function.

  The built-in method uses PLS with cross-validation (using
  \code{numSegments} random segments) to assess the
  prediction power of the variable subset. In each of the
  \code{numReplications} replications the standard error of
  prediction (SEP) is used to quantify the fitness of the
  subset. The final fitness is the mean SEP. The larger the
  number of replications, the better the estimation of the
  SEP but the slower the algorithm (the evaluation step is
  done \code{numGenerations} * \code{populationSize} times
  - see \code{\link{genAlgControl}}).

  The user specified function must take a logical vector as
  its first argument. This logical vector specifies the
  variables in the subset. The function must return a
  number representing the fitness of the variable subset
  (the higher the value the fitter the subset)
}

