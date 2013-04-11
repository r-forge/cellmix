\name{rpure}
\alias{rpure}
\title{Generating Random Pure Cell Type Sample}
\usage{
  rpure(x, n = 100, p = 20, min = 0, max = 20)
}
\arguments{
  \item{x}{number of cell types}

  \item{n}{number of features}

  \item{p}{number of samples per cell type}

  \item{min}{lower and upper limits of the distribution.
  Must be finite.}

  \item{max}{lower and upper limits of the distribution.
  Must be finite.}
}
\value{
  an \code{\link{ExpressionMix}} object that contains the
  expression data, the pure basis signatures, and a
  phenotype variable \code{'CellType'}, which indicates the
  cell type of each sample.
}
\description{
  \code{rpure} generates random expression data that
  simulates pure cell type samples, that share most gene
  expression profile pattern.
}
\examples{
x <- rpure(3)
aheatmap(x, annCol=TRUE)
}

