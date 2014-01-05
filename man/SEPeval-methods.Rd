\docType{methods}
\name{SEPeval}
\alias{SEPeval}
\alias{SEPeval,GenAlgLMEvaluator,GenAlg-method}
\alias{SEPeval,GenAlgPLSEvaluator,GenAlg-method}
\alias{SEPeval,GenAlgUserEvaluator,GenAlg-method}
\title{Get the SEP of all variable subsets}
\usage{
SEPeval(object, genAlg)

\S4method{SEPeval}{GenAlgPLSEvaluator,GenAlg}(object, genAlg)

\S4method{SEPeval}{GenAlgUserEvaluator,GenAlg}(object, genAlg)

\S4method{SEPeval}{GenAlgLMEvaluator,GenAlg}(object, genAlg)
}
\arguments{
  \item{object}{The used evaluator, i.e. an object of type
  \code{\link{GenAlgEvaluator}}}

  \item{genAlg}{An object of type \code{\link{GenAlg}}}
}
\value{
A vector with the estimated SEP for each solution
}
\description{
Calculate the standard error of prediction for all variable
subsets
}
\details{
This method is used to get the standard error of prediction
(SEP) of all variable subsets found by the genetic
algorithm.
}

