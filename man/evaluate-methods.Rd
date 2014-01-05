\docType{methods}
\name{evaluate}
\alias{evaluate}
\alias{evaluate,ANY,matrix,numeric,NULL-method}
\alias{evaluate,ANY,matrix,numeric,missing-method}
\alias{evaluate,ANY,matrix,numeric,numeric-method}
\alias{evaluate,GenAlgEvaluator,matrix,numeric,integer-method}
\alias{evaluate,GenAlgUserEvaluator,matrix,numeric,integer-method}
\title{Evaluate the fitness of variables}
\usage{
evaluate(object, X, y, seed)

\S4method{evaluate}{ANY,matrix,numeric,missing}(object, X, y, seed)

\S4method{evaluate}{ANY,matrix,numeric,NULL}(object, X, y, seed)

\S4method{evaluate}{ANY,matrix,numeric,numeric}(object, X, y, seed)

\S4method{evaluate}{GenAlgEvaluator,matrix,numeric,integer}(object, X, y, seed)

\S4method{evaluate}{GenAlgUserEvaluator,matrix,numeric,integer}(object, X, y,
  seed)
}
\arguments{
  \item{object}{The GenAlgEvaluator object that is used to
  evaluate the variables}

  \item{X}{The data matrix used to for fitting the model}

  \item{y}{The response vector}

  \item{seed}{The value to seed the random number generator
  before evaluating}
}
\description{
Evaluate the given data with the given evaluator
}

