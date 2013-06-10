\name{genAlg}
\alias{genAlg}
\title{Genetic algorithm for variable subset selection}
\usage{
  genAlg(y, X,
    control = genAlgControl(populationSize = floor(sqrt(ncol(X))), numGenerations = 100L),
    evalControl = genAlgEvalControl())
}
\arguments{
  \item{y}{The numeric response vector of length n}

  \item{X}{A n x p numeric matrix with all p covariates}

  \item{control}{Options for controlling the genetic
  algorithm. See \code{\link{genAlgControl}} for details.}

  \item{evalControl}{Options for controlling the evaluation
  step. See \code{\link{genAlgEvalControl}} for details.}
}
\value{
  Returns a list with elements \code{subsets} and
  \code{fitness}.  \item{\code{subsets}}{logical matrix
  with one variable subset per column. The columns are
  ordered according to their fitness (first column is the
  fittest)} \item{\code{fitness}}{numeric vector with the
  fitness of the corresponding variable subset}
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
  using the built-in PLS method to evaluate the standard
  error of prediction (SEP) or a user defined evaluation
  function (see \code{\link{genAlgEvalControl}}).
  Chromosomes with higher fitness have higher probability
  of mating with another chromosome. \code{populationSize /
  2} couples each create 2 children. The children are
  created by randomly mixing the parents' variables. These
  children make up the new generation and are again
  selected for mating based on their fitness. A total of
  \code{numGenerations} generations are built this way. The
  algorithm returns the last generation as well as the best
  \code{elitism} chromosomes from all generations.
}

