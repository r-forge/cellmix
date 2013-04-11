\name{.gedNNLS}
\alias{.gedNNLS}
\title{Partial Gene Expression Deconvolution by Nonnegative Least-Square}
\usage{
  .gedNNLS(X, seed, rescale = TRUE, ...)
}
\arguments{
  \item{...}{extra arguments passed to
  \code{\link[NMF]{.fcnnls}}.}

  \item{X}{target matrix}

  \item{seed}{initial NMF model as an \code{\link{NMF}}
  object.}

  \item{rescale}{logical used when estimating proportions
  from signatures, that indicates if the esti,ated
  coefficients should be scaled to sum up to one
  (\code{TRUE}) or left as estimated by the linear
  regression (\code{FALSE}). This scaling is performed
  after the coefficients have been forced to be
  nonnegative.}
}
\value{
  an \code{\link{NMF}} object.
}
\description{
  Partial Gene Expression Deconvolution by Nonnegative
  Least-Square
}
\keyword{internal}

