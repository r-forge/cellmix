\name{gedAlgorithm.nnls}
\alias{gedAlgorithm.nnls}
\alias{nnls-ged}
\title{Partial Gene Expression Deconvolution by Nonnegative Least-Squares}
\usage{
  gedAlgorithm.nnls(..., rescale = TRUE)
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
  Estimates cell/tissue proportions given a known set of
  cell/tissue-specific expression signatures, using the
  fast combinatorial nonnegative least-square method, which
  was proposed by \cite{Kim et al. (2007)} to perform
  nonnegative matrix factorization of gene expression --
  not originally for deconvolution.
}
\details{
  All regressions are fitted using the function
  \code{\link[NMF]{.fcnnls}}.
}
\examples{
# random target matrix
x <- rmatrix(100, 20)
# random cell signatures
s <- rmatrix(100, 3)

# deconvolve using nonnegative least-squares
res <- ged(x, s, 'nnls')
coef(res)
# signatures are not updated
identical(basis(res), s)
\dontshow{
	stopifnot(identical(basis(res), s))
	stopifnot( nmf.equal(res, ged(x, s, 'nnls')) )
}
}
\references{
  Kim H and Park H (2007). "Sparse non-negative matrix
  factorizations via alternating non-negativity-constrained
  least squares for microarray data analysis."
  _Bioinformatics (Oxford, England)_, *23*(12), pp.
  1495-502. ISSN 1460-2059, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btm134>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/17483501>.
}

