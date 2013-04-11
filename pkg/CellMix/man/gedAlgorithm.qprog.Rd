\name{gedAlgorithm.qprog}
\alias{cs-qprog-ged}
\alias{gedAlgorithm.cs_qprog}
\alias{gedAlgorithm.qprog}
\alias{qprog-ged}
\title{Partial Gene Expression Deconvolution by Quadratic Programming}
\usage{
  gedAlgorithm.qprog(..., data = NULL, exact = TRUE,
    log = NULL)

  gedAlgorithm.cs_qprog(..., data = NULL, exact = TRUE,
    log = NULL)
}
\arguments{
  \item{...}{extra parameters passed to
  \code{\link{mlsei}}.}

  \item{data}{markers as \code{\link{MarkerList}} object.}

  \item{exact}{logical that specifies if one should impose
  a sum-to-one (\code{TRUE}) or sum-to-less-than-one
  (\code{FALSE}) constraint on the proportions.}

  \item{log}{indicates if the data should be converted to
  linear scale (not log-trasnformed). If \code{TRUE} the
  data is exponentialised (using base 2). If \code{FALSE}
  the data is left unchanged (the detected log scale is
  displayed in verbose mode). If a number, then it is used
  as the base to exponentialise the data.}
}
\description{
  The GED method \dQuote{qprog} estimates cell/tissue
  proportions given a known set of cell/tissue-specific
  expression signatures, using quadratic programming
  techniques as proposed by \cite{Gong et al. (2011)}.

  The GED method \dQuote{cs-qprog} is an
  \strong{experimental} method that estimates
  cell/tissue-specific expression signatures given
  cell/tissue proportions, using quadratic programming
  techniques, which can optionally impose marker
  constraints on the signatures.
}
\examples{
## Example on dummy/random data
x <- rmix(3, 100, 20)
# true cell signatures
s <- basis(x)

# deconvolve using quadratic programming
res <- ged(x, s, 'qprog')
profplot(x, res)
# signatures are not updated
identical(basis(res), s)
\dontshow{
	stopifnot( identical(basis(res), s) )
	stopifnot( nmf.equal(res, ged(x, s, 'qprog')) )
}
### CS-QPROG
## Example on dummy/random data
# random target matrix
x <- rmix(3, 100, 20)
# true cell proportions
p <- coef(x)

# deconvolve using quadratic programming
res <- ged(x, p, 'cs-qprog')
profplot(basis(x), basis(res))

# signatures are not updated
identical(coef(res), p)
\dontshow{
	stopifnot( identical(coef(res), p) )
	stopifnot( nmf.equal(res, ged(x, p, 'cs-qprog')) )
}
}
\references{
  Gong T, Hartmann N, Kohane IS, Brinkmann V, Staedtler F,
  Letzkus M, Bongiovanni S and Szustakowski JD (2011).
  "Optimal deconvolution of transcriptional profiling data
  using quadratic programming with application to complex
  clinical blood samples." _PloS one_, *6*(11), pp. e27156.
  ISSN 1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0027156>, <URL:
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3217948\&tool=pmcentrez\&rendertype=abstract>.
}

