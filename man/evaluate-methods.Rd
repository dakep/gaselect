\docType{methods}
\name{evaluate}
\alias{evaluate}
\alias{evaluate,ANY,matrix,numeric,ANY,NULL-method}
\alias{evaluate,ANY,matrix,numeric,ANY,missing-method}
\alias{evaluate,ANY,matrix,numeric,ANY,numeric-method}
\alias{evaluate,ANY,matrix,numeric,logical,integer-method}
\alias{evaluate,GenAlgEvaluator,matrix,numeric,matrix,integer-method}
\title{Evaluate the fitness of variables}
\usage{
evaluate(object, X, y, subsets, seed)

\S4method{evaluate}{ANY,matrix,numeric,ANY,missing}(object, X, y, subsets, seed)

\S4method{evaluate}{ANY,matrix,numeric,ANY,NULL}(object, X, y, subsets, seed)

\S4method{evaluate}{ANY,matrix,numeric,ANY,numeric}(object, X, y, subsets, seed)

\S4method{evaluate}{ANY,matrix,numeric,logical,integer}(object, X, y, subsets,
  seed)

\S4method{evaluate}{GenAlgEvaluator,matrix,numeric,matrix,integer}(object, X, y,
  subsets, seed)
}
\arguments{
  \item{object}{The GenAlgEvaluator object that is used to
  evaluate the variables}

  \item{X}{The data matrix used to for fitting the model}

  \item{y}{The response vector}

  \item{subsets}{The logical matrix where a column stands
  for one subset to evaluate}

  \item{seed}{The value to seed the random number generator
  before evaluating}
}
\description{
Evaluate the given data with the given evaluator
}

