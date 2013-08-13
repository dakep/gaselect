\docType{methods}
\name{validData}
\alias{validData}
\alias{validData,GenAlgEvaluator,ANY-method}
\alias{validData,GenAlgLMEvaluator,ANY-method}
\alias{validData,GenAlgPLSEvaluator,ANY-method}
\title{Check if the data is valid for the evaluator}
\arguments{
  \item{object}{The evaluator}

  \item{X}{The covariates matrix}
}
\description{
  This method checks if the covariates matrix is valid for
  the evaluator
}
\examples{
evaluator <- evaluatorLM()
validData(evaluator, matrix(0, ncol = 100, nrow = 50); # not valid because more variables than observations
}

