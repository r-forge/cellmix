\name{.gedlab}
\alias{.gedlab}
\title{Gene Expression Deconvolution: Flexible Experimental Algorithm}
\usage{
  .gedlab(target, x, data = NULL, maxIter = 1000L,
    markers = c("prior+semi", "prior", "semi", "post"),
    cscale = FALSE, sscale = FALSE, ratio = 2, log = NULL,
    ...)
}
\arguments{
  \item{v}{target matrix}

  \item{x}{factorisation rank, i.e. the number of cell
  types to extract, or a complete initial NMF model.}

  \item{data}{optional marker list used to either \emph{a
  posteriori} assign estimated signatures, or \emph{a
  priori} enforce marker block patterns \strong{before}
  each iteration.}

  \item{maxIter}{maximum number of iterations}

  \item{markers}{indicates what the markers are used for:
  \describe{ \item{sQuote{prior}}{uses
  \code{\link[=DSA-ged]{DSA}} proportion estimation method
  from \cite{Zhong2013} to compute sensible initial
  proportions from average marker expression profiles in
  the mixed sample data.} \item{sQuote{semi}}{enforces
  marker block patterns after each iteration.}
  \item{sQuote{post}}{\emph{a posteriori} assigns estimated
  signatures;} }}

  \item{cscale}{logical that indicates if the estimated
  coefficient matrix, i.e. the proportions, should be
  scaled to force them to sum up to one, after each of its
  update.}

  \item{sscale}{logical that indicates if the estimated
  signature matrix, i.e. the cell-specific expression
  profiles, should be scaled to force them to have a unit
  euclidean norm, after each of its update. If
  \code{sscale=2}, then an inverse scaling is also applied
  to the coefficient matrix, so that the overall deviance
  of the model is not changed by this transformation.}

  \item{ratio}{a numeric value that indicates the maximum
  ratio allowed between marker expression values in their
  own cell type and in other cell types. It is meaningful
  only if greater than 1, but no errors is thrown if it is
  not. E.g. using \code{ratio=2} means that marker genes
  will have their expression values on signatures other
  than their own forced to a value at least twice lower
  than on their own cell types. If \code{NULL}, then marker
  expression on cell types other than their own is forced
  to zero.}

  \item{...}{extra parameters currently not used}

  \item{log}{indicates if the data are in log-scale or
  should be converted to linear-scale. This is relevant
  because the DSA algorithm assumes that the input mixed
  data are in linear scale (i.e. not log-trasnformed).

  If \code{NULL}, then data's scale is detected by
  \code{link{is_logscale}} and conversion to linear-scale
  is performed if necessary. If \code{TRUE} the data is
  exponentialised (using log base 2). If \code{FALSE} the
  data is left unchanged (the detected log scale is
  displayed in verbose mode). If a number, then it is used
  as the log base to exponentialise the data.}
}
\description{
  \code{.gedlab} is an internal function of \pkg{CellMix}
  that implements a flexible workhorse function, whose
  arguments can be tuned to generate different kind of ged
  methods. Its purpose is mainly to experiment and
  benchmark combinations of techniques, in order to design
  "better" deconvolution algorithm. For example, it is
  called by the \sQuote{deconf} algorithm (fast version).
}
\keyword{internal}

