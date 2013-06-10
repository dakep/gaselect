\name{genAlgControl}
\alias{genAlgControl}
\title{Set control arguments for the genetic algorithm}
\usage{
  genAlgControl(populationSize, numGenerations,
    minVariables, maxVariables, elitism = 10L,
    mutationProbability = 0.01, verbosity = 0L)
}
\arguments{
  \item{populationSize}{The number of "chromosomes" in the
  population (between 1 and 2^16)}

  \item{numGenerations}{The number of generations to
  produce (between 1 and 2^16)}

  \item{minVariables}{The minimum number of variables in
  the variable subset (between 0 and p - 1 where p is the
  total number of variables)}

  \item{maxVariables}{The maximum number of variables in
  the variable subset (between 1 and p, and greater than
  \code{minVariables})}

  \item{elitism}{The number of absolute best chromosomes to
  keep across all generations (between 1 and
  min(\code{populationSize} * \code{numGenerations},
  2^16))}

  \item{mutationProbability}{The probability of mutation
  (between 0 and 1)}

  \item{verbosity}{The level of verbosity. 0 means no
  output at all, 3 is very verbose.}
}
\description{
  The population must be large enough to allow the
  algorithm to explore the whole solution space. If the
  initial population is not diverse enough, the chance to
  find the global optimum is very small. Thus the more
  variables to choose from, the larger the population has
  to be.
}
\details{
  The initial population is generated randomly. Every
  chromosome uses between \code{minVariables} and
  \code{maxVariables} (uniformly distributed).

  If the mutation probability (\code{mutationProbability}
  is greater than 0, a random number of variables is
  added/removed according to a truncated geometric
  distribution to each offspring-chromosome. The resulting
  distribution of the total number of variables in the
  subset is not uniform anymore, but almost (the smaller
  the mutation probability, the more "uniform" the
  distribution). This should not be a problem for most
  applications.

  Elitism is a method of enhancing the GA by keeping track
  of very good solutions. The parameter \code{elitism}
  specifies how many "very good" solutions should be kept.
}
\seealso{
  See \code{\link{genAlgEvalControl}} for controlling the
  evaluation step during the GA
}

