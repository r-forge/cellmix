\name{.gedLSfit}
\alias{.gedLSfit}
\title{Partial Gene Expression Deconvolution by Least-Square}
\usage{
  .gedLSfit(X, seed, nneg = c("iterate", "pmax", "none"),
    rescale = TRUE, ...)
}
\arguments{
  \item{nneg}{specification of the method used to enforce
  the nonnegativity of the estimated proportions (used only
  in method \dQuote{lsfit}).}

  \item{rescale}{logical used when estimating proportions
  from signatures, that indicates if the esti,ated
  coefficients should be scaled to sum up to one
  (\code{TRUE}) or left as estimated by the linear
  regression (\code{FALSE}). This scaling is performed
  after the coefficients have been forced to be
  nonnegative.}

  \item{...}{extra arguments passed to
  \code{\link{lsfit}}.}

  \item{X}{target matrix}

  \item{seed}{initial NMF model as an \code{\link{NMF}}
  object.}
}
\value{
  an \code{\link{NMF}} object.
}
\description{
  Partial Gene Expression Deconvolution by Least-Square
}
\keyword{internal}

