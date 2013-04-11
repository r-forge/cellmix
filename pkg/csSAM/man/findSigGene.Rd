\name{findSigGene}
\alias{findSigGene}
\title{findSigGene}
\usage{
  findSigGene(G, cc, y, rhat, csSAMData)
}
\arguments{
  \item{G}{Gene expression matrix of heterogenous tissue
  measurements}

  \item{cc}{Matrix of cell-frequency measures per person}

  \item{y}{Numeric group association of each sample. Either
  1 or 2.}

  \item{rhat}{Matrix of cell-specific contrasts for each
  gene in each cell-type as computed for the original group
  classification.}

  \item{csSAMData}{List object returned from fdrCsSAM.}
}
\value{
  A matrix size k by g where k is the number of cell-types
  and g is the number of genes. For each cell in the
  matirx, listed is the FDR of the gene for a difference in
  a given cell-type.
}
\description{
  Find the false discovery rate for each gene in each
  cell-type.
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

