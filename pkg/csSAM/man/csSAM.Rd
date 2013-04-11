\name{csSAM}
\alias{csSAM}
\title{csSAM}
\usage{
  csSAM(ghat1, se1, n1, ghat2, se2, n2, standardize,
    medianCenter = TRUE, nonNeg = FALSE)
}
\arguments{
  \item{ghat1}{Expression matrix of deconvolved
  cell-specific gene expression estimates for group 1.}

  \item{se1}{Standard error group 1}

  \item{n1}{Group 1 size}

  \item{ghat2}{Expression matrix of deconvolved
  cell-specific gene expression estimates for group 2.}

  \item{se2}{Standard error group 2}

  \item{n2}{Group 2 size}

  \item{standardize}{Standardize contrast values}

  \item{medianCenter}{Median center rhat distributions for
  each cell-type}

  \item{nonNeg}{Negative values not allowed such as in a
  single channel microarray. Zero them if negative (a
  conervative option)}
}
\value{
  A matrix object with the result of contrasting the
  average cell-specific expression profile of the two
  groups, per cell-type (Size k by g where k is the number
  of cells and g is the number of genes).
}
\description{
  Computes the constrast between groups for the deconvolved
  cell-specific expression for each cell-type
}
\author{
  Shai Shen-Orr, Rob Tibshirani, Narasimhan
  Balasubramanian, David Wang
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

