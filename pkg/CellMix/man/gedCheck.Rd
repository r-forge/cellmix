\name{gedCheck}
\alias{gedCheck}
\title{Checking a Deconvolution Method}
\usage{
  gedCheck(method, ..., maxIter = 5L, verbose = 3L)
}
\arguments{
  \item{method}{algorithm specification}

  \item{...}{other arguments passed to \code{\link{ged}}.}

  \item{maxIter}{maximum number of iterations. Note that it
  is fixed to a small default value for the purpose of the
  check, which is likely to make no sense for real
  application.}

  \item{verbose}{defines verbosity level (\code{FALSE/TRUE}
  or a number).}
}
\value{
  returns invisibly the result of running the algorithm,
  which usually is an \code{NMFfit} object.
}
\description{
  This function runs a given algorithm on a small random
  toy dataset, to check that the algorithm works correctly
}
\examples{
gedCheck('deconf')
gedCheck('DSA', log=FALSE)
}

