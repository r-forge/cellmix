\name{plotCsSAM}
\alias{plotCsSAM}
\title{plotcsSAM}
\usage{
  plotCsSAM(csSAMdata, SAMdata, alternative, cellID,
    numcell, fileName)
}
\arguments{
  \item{csSAMdata}{List object output of the fdrCsSAM
  function}

  \item{SAMdata}{List object output of the fdrSAM function}

  \item{alternative}{Type of test conducted. Will appear in
  plot title.}

  \item{cellID}{Label for each cell-type}

  \item{numcell}{Number of different cell-types being
  considered.}

  \item{fileName}{Name of output pdf file.}
}
\description{
  Plots the # of genes called significnat at a given false
  disocvery rate for the SAM (heterogenous tissue)
  comparison, and for each of the contrasted cell-types
  using csSAM
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

