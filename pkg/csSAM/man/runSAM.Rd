\name{runSAM}
\alias{runSAM}
\title{runSAM}
\usage{
  runSAM(G, y, s0.sam = NULL, stand.r = TRUE)
}
\arguments{
  \item{G}{Matrix of gene expression, columns ordered in
  the same order at the cell-frequency matrix (n by p, n
  samples, p genes)}

  \item{y}{Numeric group association of each sample. Either
  1 or 2.}

  \item{s0.sam}{Input or computed value of SAM
  exchangeability factor. Default is determined
  automatically}

  \item{stand.r}{Median center and standardize arrays.
  Default is TRUE.}
}
\description{
  A lightweight version of the SAM algorithm, only performs
  two group comparison with equal deltas on each tail
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

