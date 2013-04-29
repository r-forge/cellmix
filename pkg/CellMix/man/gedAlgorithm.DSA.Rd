\name{gedAlgorithm.DSA}
\alias{DSA-ged}
\alias{gedAlgorithm.DSA}
\title{Complete Deconvolution using Digital Sorting Algorithm (DSA)}
\usage{
  gedAlgorithm.DSA(..., data = NULL, exact = TRUE,
    log = NULL)
}
\arguments{
  \item{...}{extra parameters DSAproportionspassed to
  \code{\link{mlsei}}.}

  \item{data}{markers as \code{\link{MarkerList}} object.}

  \item{exact}{logical that specifies if one should impose
  a sum-to-one (\code{TRUE}) or sum-to-less-than-one
  (\code{FALSE}) constraint on the proportions.}

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
  The method \dQuote{DSA} implements the Digital Sorting
  Algorithm (DSA) proposed by \cite{Zhong et al. (2013)},
  which performs complete gene expression deconvolution
  using a set of marker genes only.
}
\examples{
## Example on dummy/random data
x <- rmix(3, 100, 20)
# markers
ml <- getMarkers(x)

# deconvolve using DSA (quadratic programming)
res <- ged(x, ml, 'DSA', verbose=TRUE, log=FALSE)
# plot against true data
profplot(x, res)
profplot(t(basis(x)), t(basis(res)), legend=FALSE)

\dontshow{
	stopifnot( nmf.equal(res, ged(x, ml, 'DSA', log=FALSE), identical=FALSE) )
}
}
\references{
  Zhong Y, Wan Y, Pang K, Chow LM and Liu Z (2013).
  "Digital sorting of complex tissues for cell
  type-specific gene expression profiles." _BMC
  Bioinformatics_, *14*(1), pp. 89. ISSN 1471-2105, <URL:
  http://dx.doi.org/10.1186/1471-2105-14-89>, <URL:
  http://www.biomedcentral.com/1471-2105/14/89>.
}

