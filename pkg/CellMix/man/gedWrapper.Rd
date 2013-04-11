\name{gedWrapper}
\alias{gedWrapper}
\title{Wrapping GED Algorithms}
\usage{
  gedWrapper(key)
}
\arguments{
  \item{key}{Name of the GED algorithm to be wrapped. It
  should be the name of a registered algorithm as returned
  by \code{\link{gedAlgorithm}}, or a GED algorithm object
  (i.e. an instance of \code{\linkS4class{GEDStrategy}}).}
}
\value{
  a function with attribute \code{'algorithm'} set to the
  value of \code{key}
}
\description{
  This function creates a wrapper function for calling the
  function \code{\link{ged}} with a given GED algorithm.
}
\examples{
# wrap 'qprog' algorithm from Gong et al. (2011)
qp <- nmfWrapper('qprog')
qp

# test on random data
x <- rmatrix(100,20)
sig <- rmatrix(100,3)
res <- ged(x, sig, 'qprog', seed=12345)
res2 <- qp(x, sig, seed=12345)
nmf.equal(res, res2)

\dontshow{ stopifnot(nmf.equal(res, res2)) }
}
\seealso{
  \code{\link{gedAlgorithm}}, \code{\link{ged}}
}
\keyword{internal}

