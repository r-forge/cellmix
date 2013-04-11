\name{gedAlgorithm.deconf}
\alias{.deconf-ged}
\alias{deconf-ged}
\alias{gedAlgorithm.deconf}
\alias{gedAlgorithm._deconf}
\title{Complete Gene Expression Deconvolution: Method deconf}
\usage{
  gedAlgorithm.deconf(..., maxIter = 1000L, data = NULL,
    error.threshold = 0)

  gedAlgorithm._deconf(..., maxIter = 1000L, data = NULL,
    error.threshold = 0)
}
\arguments{
  \item{...}{extra parameters currently not used}

  \item{maxIter}{maximum number of iterations}

  \item{data}{optional marker list used to either \emph{a
  posteriori} assign estimated signatures, or \emph{a
  priori} enforce marker block patterns \strong{before}
  each iteration.}

  \item{error.threshold}{ min error to continue iterations,
  if error falls below this threshold, the iterative
  algorithm stops }
}
\description{
  The \code{\link{ged}} method \sQuote{deconf} uses an
  alternate least-squares algorithm to estimate both cell
  proportions and cell-specific signatures from global
  expression data, as proposed by \cite{Repsilber et al.
  (2010)}.
}
\details{
  This method fits an NMF model to the data in a completely
  \emph{unsupervised} manner. If marker genes are provided,
  they are used \strong{a posteriori} to assign each
  estimated component, i.e. each cell-specific signature,
  to the cell-type with the greatest proportions of
  consistent markers.
}
\section{Fast built-in implementation}{
  The method \sQuote{deconf} is implemented as an NMF
  algorithm, which is registered under the same names in
  the \pkg{NMF} package's algorithm registry.

  It uses an improved implementation, based on the fast
  combinatorial nonnegative least-squares algorithm from
  \cite{Van Benthem et al. (2004), Kim et al. (2007)}, as provided by the
  function \code{\link{fcnnls}} in the \pkg{NMF} package.
  This enables to achieve great performance speed-up, being
  really -- way -- much faster than the original
  implementation.
}

\section{Original implementation}{
  The \pkg{CellMix} also includes a way to run the original
  version from the \pkg{deconf} package, registered under
  the name \sQuote{.deconf} in the algorithm registry.

  This version requires the \pkg{deconf} package, which was
  released as supplementary data only to support the paper
  from \cite{Repsilber et al. (2010)}, i.e. it is not available from
  CRAN or Bioconductor. However, we made it available from
  the \pkg{CellMix} CRAN-like support repository:

  \url{http://web.cbio.uct.ac.za/~renaud/CRAN}

  The easiest way to install it is to run:

  \code{install.extras('CellMix', 'deconf')}
}
\references{
  Repsilber D, Kern S, Telaar A, Walzl G, Black GF, Selbig
  J, Parida SK, Kaufmann SHE and Jacobsen M (2010).
  "Biomarker discovery in heterogeneous tissue samples
  -taking the in-silico deconfounding approach." _BMC
  bioinformatics_, *11*, pp. 27. ISSN 1471-2105, <URL:
  http://dx.doi.org/10.1186/1471-2105-11-27>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/20070912>.

  Van Benthem M and Keenan MR (2004). "Fast algorithm for
  the solution of large-scale non-negativity-constrained
  least squares problems." _Journal of Chemometrics_,
  *18*(10), pp. 441-450. ISSN 0886-9383, <URL:
  http://dx.doi.org/10.1002/cem.889>, <URL:
  http://doi.wiley.com/10.1002/cem.889>.

  Kim H and Park H (2007). "Sparse non-negative matrix
  factorizations via alternating non-negativity-constrained
  least squares for microarray data analysis."
  _Bioinformatics (Oxford, England)_, *23*(12), pp.
  1495-502. ISSN 1460-2059, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btm134>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/17483501>.
}

