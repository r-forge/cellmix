\name{gedAlgorithm.ssKL}
\alias{gedAlgorithm.ssFrobenius}
\alias{gedAlgorithm.ssKL}
\alias{ssFrobenius-ged}
\alias{ssKL-ged}
\title{Complete Gene Expression Deconvolution by Semi-Supervised NMF}
\usage{
  gedAlgorithm.ssKL(..., maxIter = 3000L, seed = "rprop",
    eps = 2.22044604925031e-16, .stop = NULL, data = NULL,
    markers = c("prior+semi", "semi"), ratio = NULL,
    copy = FALSE, sscale = FALSE, alpha = 0,
    stationary.th = .Machine$double.eps,
    check.interval = 5 * check.niter, check.niter = 10L)

  gedAlgorithm.ssFrobenius(..., maxIter = 3000L,
    seed = "rprop", .stop = NULL, data = NULL,
    markers = c("prior+semi", "semi"), ratio = NULL,
    eps = NULL, sscale = TRUE, copy = FALSE, alpha = 0,
    stationary.th = .Machine$double.eps,
    check.interval = 5 * check.niter, check.niter = 10L)
}
\arguments{
  \item{seed}{default seeding method.}

  \item{sscale}{specifies how signatures -- and proportions
  -- are re-scaled at the end of each iteration. If
  \code{TRUE}, each signature is mean-centered separately.
  If \code{2}, then each signature is mean-centered
  separately and the inverse linear transformation
  proportions is applied to the proportions (i.e. on the
  rows of the mixture coefficient matrix), so that the
  fitted matrix does not change. If \code{FALSE}, no
  re-scaling is performed at all.}

  \item{alpha}{numeric coefficient used to smoothly enforce
  a sum-up-to-one constraint on the proportions, by
  regularising the objective function. If \code{NULL,} no
  constraint is applied.}

  \item{...}{extra arguments passed to the function
  \code{\link{objective}}, which computes the objective
  value between \code{x} and \code{y}.}

  \item{stationary.th}{maximum absolute value of the
  gradient, for the objective function to be considered
  stationary.}

  \item{check.interval}{interval (in number of iterations)
  on which the stopping criterion is computed.}

  \item{check.niter}{number of successive iteration used to
  compute the stationnary criterion.}

  \item{maxIter}{maximum number of iterations to perform.}

  \item{eps}{small numeric value used to ensure numeric
  stability, by shifting up entries from zero to this fixed
  value.}

  \item{.stop}{specification of a stopping criterion, that
  is used instead of the one associated to the NMF
  algorithm. It may be specified as: \itemize{ \item the
  access key of a registered stopping criterion; \item a
  single integer that specifies the exact number of
  iterations to perform, which will be honoured unless a
  lower value is explicitly passed in argument
  \code{maxIter}. \item a single numeric value that
  specifies the stationnarity threshold for the objective
  function, used in with \code{\link{nmf.stop.stationary}};
  \item a function with signature
  \code{(object="NMFStrategy", i="integer", y="matrix",
  x="NMF", ...)}, where \code{object} is the
  \code{NMFStrategy} object that describes the algorithm
  being run, \code{i} is the current iteration, \code{y} is
  the target matrix and \code{x} is the current value of
  the NMF model. }}

  \item{copy}{logical that indicates if the update should
  be made on the original matrix directly (\code{FALSE}) or
  on a copy (\code{TRUE} - default). With \code{copy=FALSE}
  the memory footprint is very small, and some speed-up may
  be achieved in the case of big matrices. However, greater
  care should be taken due the side effect. We recommend
  that only experienced users use \code{copy=TRUE}.}

  \item{data}{marker list}

  \item{markers}{indicates what the markers are used for:
  \describe{ \item{sQuote{prior}}{uses
  \code{\link[=DSA-ged]{DSA}} proportion estimation method
  from \cite{Zhong2013} to compute sensible initial
  proportions from average marker expression profiles in
  the mixed sample data.} \item{sQuote{semi}}{enforces
  marker block patterns after each iteration.}
  \item{sQuote{post}}{\emph{a posteriori} assigns estimated
  signatures;} }}

  \item{ratio}{expression ratio of markers between its cell
  type and other cell types.}
}
\description{
  Algorithms \sQuote{ssKL} and \sQuote{ssFrobenius} are
  modified versions of the original NMF algorithm from
  \cite{Brunet et al. (2004)} and \cite{Lee et al. (2001)},
  that use a set of known marker genes for each cell type,
  to enforce the expected block expression pattern on the
  estimated signatures, as proposed in \cite{Gaujoux et al.
  (2011)}.
}
\details{
  These algorithms simultaneously estimates both the
  cell-specific signature and mixture proportion matrices,
  using block-descent method that alternately estimates
  each matrix. Both re-scale the final proportion estimates
  so that they sum-up to one.

  The functions \code{gedAlgorithm.ssKL} and
  \code{gedAlgorithm.ssFrobenius} are wrapper functions to
  the underlying NMF algorithms. They are primiraly defined
  to enable correct listing their specific arguments on
  this help page. The recommend way of applying these
  algorithms is via \code{\link{ged}} interface (e.g.,
  \code{ged(..., method='ssKL')}).
}
\references{
  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.

  Lee DD and Seung H (2001). "Algorithms for non-negative
  matrix factorization." _Advances in neural information
  processing systems_. <URL:
  http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.

  Gaujoux R and Seoighe C (2011). "Semi-supervised
  Nonnegative Matrix Factorization for gene expression
  deconvolution: A case study." _Infection, genetics and
  evolution : journal of molecular epidemiology and
  evolutionary genetics in infectious diseases_. ISSN
  1567-7257, <URL:
  http://dx.doi.org/10.1016/j.meegid.2011.08.014>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/21930246>.
}

