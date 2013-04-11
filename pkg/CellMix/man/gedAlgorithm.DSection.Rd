\name{gedAlgorithm.DSection}
\alias{csTopTable.ged_DSection}
\alias{DSection-ged}
\alias{gedAlgorithm.DSection}
\title{Partial Gene Expression Deconvolution with DSection}
\usage{
  gedAlgorithm.DSection(..., maxIter = 500L, p0,
    data = NULL, nBurnIn = 4 * maxIter,
    W0 = max(100, nbasis(seed)), W_proposal = 100,
    samplep = TRUE)

  \method{csTopTable}{ged_DSection} (x, coef = 1L,
    decreasing = TRUE, ...)
}
\arguments{
  \item{data}{variable (e.g., factor) that defines groups
  of samples. Cell-specific signatures will be computed for
  each group.}

  \item{maxIter}{number of sampling, not including burn-in
  sampling.}

  \item{...}{extra arguments passed to
  \code{\link{DSection}}.}

  \item{coef}{specifies the reference cell type or group of
  samples, for which differential expression is computed.}

  \item{p0}{<T-by-J> matrix of prior predictions on cell
  type proportions. T is the number of cell types, and
  columns in \code{p0} must be positive and add up to one.}

  \item{nBurnIn}{Amount of burn-in. NOTE: keep nBurnIn >
  0.}

  \item{W0}{Prior prediction weight, i.e., degree of
  confidence, on \code{p0}. Defines the peakedness of
  Dirichlet density around p0. NOTE: keep W0 >= T.}

  \item{W_proposal}{Transition kernel weight, defines the
  peakedness of Dirichlet density around p*, the old value.
  The higher W_proposal is, the smaller the proposal steps
  around p* are.}

  \item{samplep}{logical value, indicating whether to
  sample from the posterior for cell type proportions
  (\code{TRUE}) or not (\code{FALSE}). SUGGESTED USE:
  sample from the posterior (samplep = 1).}

  \item{x}{data object, typically returned by
  \code{\link{ged}}.}

  \item{decreasing}{logical that indicates the feature
  ordering, based on their p-values or FDRs.}
}
\description{
  Estimates cell/tissue cell/tissue-specific expression
  signatures, given proportion priors using the MCMC
  approach from \cite{Erkkila et al. (2010)}, implemented
  in Matlab and wrapped in \pkg{CellMix} by the function
  \code{\link{DSection}}.

  The S3 method \code{csTopTable} for DSection fits
  computes nominal p-values (i.e. unadjusted) of
  differential expression between cell type or group of
  samples within each cell type, for deconvolution results
  from the \code{\link[=DSection-ged]{DSection}} algorithm.
}
\examples{
# random global expression
x <- rmix(3, 50, 20)
# extract true proportions
p <- coef(x)
# add noise
p0 <- scoef(abs(p + rmatrix(p, dist=rnorm, sd=0.15)))
# check how noisy this makes the proportion estimates
profplot(p, p0)

# deconvolve using DSection
\dontrun{
res <- ged(x, p0, 'DSection', maxIter=10, seed=12345)
head(basis(res))
# proportions are updated
!identical(coef(res), p0)
# check how better they are
profplot(x, res)

\dontshow{
	stopifnot( !identical(coef(res), p0) )
	stopifnot( nmf.equal(res, ged(x, p0, 'DSection', maxIter=10, seed=12345)) )
}
}
}
\references{
  Erkkila T, Lehmusvaara S, Ruusuvuori P, Visakorpi T,
  Shmulevich I and Lahdesmaki H (2010). "Probabilistic
  analysis of gene expression measurements from
  heterogeneous tissues." _Bioinformatics (Oxford,
  England)_, *26*(20), pp. 2571-7. ISSN 1367-4811, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btq406>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/20631160
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2951082\&tool=pmcentrez\&rendertype=abstract>.
}

