\docType{methods}
\name{log,ExpressionSet-method}
\alias{expb}
\alias{expb,ExpressionSet-method}
\alias{expb,matrix-method}
\alias{expb-methods}
\alias{exp,ExpressionSet-method}
\alias{log,ExpressionSet-method}
\alias{quantile.ExpressionSet}
\alias{range,ExpressionSet-method}
\title{Numeric Computations on ExpressionSet objects}
\usage{
  \S4method{log}{ExpressionSet}(x, ...)

  expb(x, ...)

  \S4method{expb}{matrix}(x, base = exp(1))

  \S4method{exp}{ExpressionSet}(x)

  \S4method{range}{ExpressionSet}(x, ..., na.rm = FALSE)

  \method{quantile}{ExpressionSet} (x, ...)
}
\arguments{
  \item{x}{an \code{ExpressionSet} object.}

  \item{...}{extra arguments passed to subsequent calls,
  usually of the corresponding method in the \pkg{stats}
  package.}

  \item{base}{log base to use.}

  \item{na.rm}{logical, indicating if \code{\link{NA}}'s
  should be omitted.}
}
\description{
  The \pkg{CellMix} defines some generics and methods to
  apply numeric transformations to \code{ExpressionSet}
  objects, which is convenient when working on gene
  expression deconvolution algorithms, where scale
  (log/linear) may matter.

  \code{log} log-transforms the expression matrix of
  \code{\link{ExpressionSet}} objects.

  \code{expb} applies an entry-wise exponential
  transformation to the expression matrix of
  \code{\link{ExpressionSet}}, in a specific base.

  \code{quantile} computes the range of expression values
  in \code{\link{ExpressionSet}} objects.
}
\section{Methods}{
  \describe{

  \item{exp}{\code{signature(x = "ExpressionSet")}:
  \code{exp} is equivalent to \code{expb(x, 1)}. }

  \item{range}{\code{signature(x = "ExpressionSet")}:
  \code{range} computes the range of expression values from
  \code{\link{ExpressionSet}} objects. }

  }
}
\keyword{methods}

