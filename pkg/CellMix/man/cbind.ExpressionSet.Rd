\name{cbind.ExpressionSet}
\alias{cbind.ExpressionSet}
\title{Combining Expression Matrices}
\usage{
  \method{cbind}{ExpressionSet} (..., deparse.level = 1)
}
\arguments{
  \item{...}{series of \code{ExpressionSet} and/or matrix
  objects.}

  \item{deparse.level}{integer controlling the construction
  of labels in the case of non-matrix-like arguments (for
  the default method):\cr \code{deparse.level = 0}
  constructs no labels; the default,\cr \code{deparse.level
  = 1 or 2} constructs labels from the argument names, see
  the \sQuote{Value} section below.}
}
\value{
  an \code{ExpressionSet} object
}
\description{
  The method \code{cbind.ExpressionSet} combines sample
  expression profiles from multiple
  \code{\link{ExpressionSet}} or matrix objects.
}
\details{
  The expression matrices must be exactly of the same
  dimensions. For the end result to be meaningful, one
  probably wants the row names to match as well, i.e. both
  object contain the same features, in the same order.
  However no check is done for this.

  Note that the returned \code{ExpressionSet} object has no
  sample or feature annotations.
}

