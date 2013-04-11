\name{gedAlgorithm.csSAM}
\alias{csplot.ged_csSAM}
\alias{csSAM-ged}
\alias{csTopTable.ged_csSAM}
\alias{gedAlgorithm.csSAM}
\title{Partial Gene Expression Deconvolution with csSAM}
\usage{
  gedAlgorithm.csSAM(Y, x, data = NULL, nperms = 200,
    alternative = c("two.sided", "greater", "less"),
    standardize = TRUE, medianCenter = TRUE, logRm = FALSE,
    logBase = 2, nonNeg = TRUE, verbose = lverbose())

  \method{csTopTable}{ged_csSAM} (x, ...)

  \method{csplot}{ged_csSAM} (x, types = NULL,
    xlab = "# called", ylab = "FDR", ylim = c(0, 1), ...)
}
\arguments{
  \item{types}{index or names of the type to plot. They
  need to be found in the fit data.}

  \item{Y}{target global gene expression matrix (n x p),
  with samples in columns, ordered in the same order at the
  cell proportions data in \var{x}.}

  \item{x}{known cell proportions as a matrix (k x p) or an
  \code{\linkS4class{NMF}} model containing the cell
  proportions in the coefficient matrix -- and a normally
  empty basis matrix. The proportions must be ordered in
  the same order as the samples in the target matrix.

  For \code{csTopTable}, a csSAM fit as return by
  \code{\link{ged}}.}

  \item{data}{specification of the sample groups. If not
  missing, it must be a factor or coercible to a factor,
  with length the number of samples, i.e. columns, in the
  target matrix.}

  \item{nperms}{The number of permutations to perform. It
  is only used when computing cell-specific differential
  expression between groups specified in argument
  \code{data}.}

  \item{alternative}{two.sided less greater}

  \item{standardize}{Standardize sample or not. Default is
  TRUE.}

  \item{medianCenter}{Median center rhat distributions.
  Default is TRUE.}

  \item{logRm}{Exponentiate data for deconvolution stage.
  Default is FALSE}

  \item{logBase}{Base of logaritm used to determine
  exponentiation factor. Default is 2}

  \item{nonNeg}{For single channel arrays. Set any
  cell-specific expression estimated as negative, to a
  ceiling of 0. It is conservative in its study of
  differential expression. Default is FALSE.}

  \item{verbose}{logical that indicates if verbose messages
  should be shown.}

  \item{...}{extra parameters passed to subsequent calls.}

  \item{xlab}{a label for the x axis, defaults to a
  description of \code{x}.}

  \item{ylab}{a label for the y axis, defaults to a
  description of \code{y}.}

  \item{ylim}{the y limits of the plot.}
}
\description{
  Estimates cell/tissue proportions given a known set of
  cell/tissue-specific expression signatures, using
  standard least-squares, as implemented by the package
  \code{\link{csSAM}}.

  The S3 method \code{csTopTable} for csSAM fits returns,
  for each feature, the false discovery rates of
  differential expression between groups of samples within
  each cell type, as computed by
  \code{\link[=csSAM]{fdrCsSAM}} when running csSAM. These
  are returned as a list, whith one element per cell type.

  The S3 method \code{csplot} for csSAM fits plots
  cell-specific fdr cumulative distributions.
}
\details{
  All regressions are fitted using the function
  \code{\link{lsfit}}.
}
\examples{
# random global expression
x <- rmix(3, 100, 20)
basisnames(x) <- paste('Cell', 1:nbasis(x))
# extract true proportions
p <- coef(x)

# deconvolve using csSAM
res <- ged(x, p, 'csSAM')
head(basis(res))
# proportions are not updated
identical(coef(res), p)
\dontshow{
	stopifnot(identical(coef(res), p))
	stopifnot( nmf.equal(res, ged(x, p, 'csSAM')) )
}

# estimate cell-specific differential expression between 2 groups
gr <- gl(2, 10)
res <- ged(x, p, 'csSAM', data = gr, nperms=20, verbose=TRUE)
head(basis(res))
# plot FDRs
csplot(res)
# extract fdr for top differentially expressed gene in each cell type
t <- csTopTable(res)
str(t)
}
\references{
  Shen-Orr SS, Tibshirani R, Khatri P, Bodian DL, Staedtler
  F, Perry NM, Hastie T, Sarwal MM, Davis MM and Butte AJ
  (2010). "Cell type-specific gene expression differences
  in complex tissues." _Nature methods_, *7*(4), pp. 287-9.
  ISSN 1548-7105, <URL:
  http://dx.doi.org/10.1038/nmeth.1439>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/20208531>.
}
\seealso{
  \code{\link[csSAM]{fdrCsSAM}}, \code{\link{csTopTable}}
}

