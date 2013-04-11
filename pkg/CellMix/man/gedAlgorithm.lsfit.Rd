\name{gedAlgorithm.lsfit}
\alias{gedAlgorithm.lsfit}
\alias{lsfit-ged}
\title{Partial Gene Expression Deconvolution by Standard Least-Squares}
\usage{
  gedAlgorithm.lsfit(...,
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
  Estimates cell/tissue proportions given a known set of
  cell/tissue-specific expression signatures, using
  standard least-squares as proposed by \cite{Abbas et al.
  (2009)}.
}
\details{
  The algorithm uses a heuristic to enforce the
  nonnegativity of the estimated proportions, that consists
  in fitting successive regressions, each time excluding
  the most negative coefficient from the model, until all
  coefficients are nonnegative.

  All regressions are fitted using the function
  \code{\link{lsfit}}.
}
\examples{
# random target matrix
x <- rmatrix(100, 20)
# random cell signatures
s <- rmatrix(100, 3)

# deconvolve using standard least-squares
res <- ged(x, s, 'lsfit')
coef(res)
# signatures are not updated
identical(basis(res), s)
\dontshow{
	stopifnot(identical(basis(res), s))
	stopifnot( nmf.equal(res, ged(x, s, 'lsfit')) )
}
}
\references{
  Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z and Clark
  HF (2009). "Deconvolution of blood microarray data
  identifies cellular activation patterns in systemic lupus
  erythematosus." _PloS one_, *4*(7), pp. e6098. ISSN
  1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0006098>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/19568420>.
}

