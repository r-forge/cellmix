\name{gedAlgorithm.cs_nnls}
\alias{cs-nnls-ged}
\alias{gedAlgorithm.cs_nnls}
\title{Partial Cell-Specific Expression by Nonnegative Least-Squares}
\usage{
  gedAlgorithm.cs_nnls(..., rescale = TRUE)
}
\arguments{
  \item{...}{extra arguments passed to
  \code{\link[NMF]{.fcnnls}}.}

  \item{rescale}{logical used when estimating proportions
  from signatures, that indicates if the esti,ated
  coefficients should be scaled to sum up to one
  (\code{TRUE}) or left as estimated by the linear
  regression (\code{FALSE}). This scaling is performed
  after the coefficients have been forced to be
  nonnegative.}
}
\description{
  Estimates cell-specific signatures given known cell
  proportions, using the fast combinatorial nonnegative
  least-square as the ged algorithm
  \code{\link[=nnls-ged]{nnls}}. It is included in the
  \pkg{CellMix} package for test/experimental purposes.
}
\details{
  All regressions are fitted using the function
  \code{\link[NMF]{.fcnnls}}.
}
\examples{
# random target matrix
x <- rmatrix(100, 20)
# random cell proprtions
p <- rmatrix(3, 20)

# deconvolve using nonnegative least-squares
res <- ged(x, p, 'cs-nnls')
head(basis(res))
# proportions are not updated
identical(coef(res), p)
\dontshow{
	stopifnot(identical(coef(res), p))
	stopifnot( nmf.equal(res, ged(x, p, 'cs-nnls')) )
}
}

