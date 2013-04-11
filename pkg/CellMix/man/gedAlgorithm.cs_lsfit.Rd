\name{gedAlgorithm.cs_lsfit}
\alias{cs-lsfit-ged}
\alias{gedAlgorithm.cs_lsfit}
\title{Cell-Specific Expression by Standard Least-Squares}
\usage{
  gedAlgorithm.cs_lsfit(...,
    nneg = c("iterate", "pmax", "none"), rescale = TRUE)
}
\arguments{
  \item{...}{extra arguments passed to
  \code{\link{lsfit}}.}

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
}
\description{
  Estimates cell-specific proportions given known
  proportions expression signatures, using standard
  least-squares.
}
\details{
  The algorithm uses the same heuristic as the ged
  algorithm \code{\link[=lsfit-ged]{lsfit}} to enforce
  non-negative expression values. It is included in the
  \pkg{CellMix} package for test/experimental purposes.
}
\examples{
# random target matrix
x <- rmatrix(100, 20)
# random cell proprtions
p <- rmatrix(3, 20)

# deconvolve using standard least-squares
res <- ged(x, p, 'cs-lsfit')
head(basis(res))
# proportions are not updated
identical(coef(res), p)
\dontshow{
	stopifnot(identical(coef(res), p))
	stopifnot( nmf.equal(res, ged(x, p, 'cs-lsfit')) )
}
}

