\docType{methods}
\name{getEvalFun}
\alias{getEvalFun}
\alias{getEvalFun,GenAlgEvaluator,GenAlg-method}
\alias{getEvalFun,GenAlgUserEvaluator,GenAlg-method}
\title{Get the evaluation function from a GenAlgUserEvaluator}
\arguments{
  \item{object}{The evaluator}

  \item{genAlg}{The \code{\link{GenAlg}} object}
}
\description{
  This method returns the correct evaluation function from
  a GenAlgUserEvaluator that can be used by the C++-code as
  callback or NULL for any other evaluator
}

