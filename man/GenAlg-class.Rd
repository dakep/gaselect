\name{GenAlg}
\alias{GenAlg}
\title{Result of a genetic algorithm run}
\description{
  Return object of a run of the genetic algorithm genAlg
}
\section{Slots}{
  \describe{ \item{\code{subsets}:}{Logical matrix with one
  variable subset per column. The columns are ordered
  according to their fitness (first column contains the
  fittest variable-subset).}
  \item{\code{rawFitness}:}{Numeric vector with the raw
  fitness of the corresponding variable subset returned by
  the evaluator.} \item{\code{response}:}{The original
  response vector.} \item{\code{covariates}:}{The original
  covariates matrix.} \item{\code{evaluator}:}{The
  evaluator used in the genetic algorithm.}
  \item{\code{control}:}{The control object.} }
}

