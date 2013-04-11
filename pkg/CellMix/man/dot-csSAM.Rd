\name{.csSAM}
\alias{.csSAM}
\title{Cell-specific Differential Expression with csSAM}
\usage{
  .csSAM(Y, x, data = NULL, nperms = 200,
    alternative = c("two.sided", "greater", "less"),
    standardize = TRUE, medianCenter = TRUE, logRm = FALSE,
    logBase = 2, nonNeg = TRUE, verbose = lverbose())
}
\arguments{
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

  \item{verbose}{logical that indicates if verbose messages
  should be shown.}

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
}
\value{
  Returns an NMF model object, with the cell specific
  signatures stored in the \code{\link{basis}} matrix, and
  the following miscellaneous slots that store details
  about the csSAM fit (accessible via \code{$}):
  \item{deconv}{A list object containing a fit (cell-type
  specific expression) for each group. Each element in the
  list is an object returned by csFit.} \item{fdr.csSAM}{ A
  list output of the fdrCsSAM function.}
  \item{sigGene.csSAM}{A list of significant genes.}
}
\description{
  This function is adapted from the function
  \code{\link[csSAM]{csSamWrapper}} in the
  \code{\link{csSAM}} package, to integrate the csSAM
  algorithm from \cite{Shen-Orr et al. (2010)} into the
  \pkg{CellMix} framework of deconvolution algorithms.
}
\author{
  Original function: Shai Shen-Orr, Rob Tibshirani,
  Narasimhan Balasubramanian, David Wang

  Adaptation for \pkg{CellMix}: Renaud Gaujoux
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
  \code{\link{csfit}},\code{\link{csSAM}},\code{\link{fdrCsSAM}}
}

