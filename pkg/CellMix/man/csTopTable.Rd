\name{csTopTable}
\alias{csTopTable}
\alias{csTopTable.array}
\alias{csTopTable.matrix}
\alias{csTopTable.NMFfit}
\alias{csTopTable.NMFfitX}
\title{Compute Cell-Specific Statistics}
\usage{
  csTopTable(x, ...)

  \method{csTopTable}{matrix} (x, n = 100L, ct = NULL,
    decreasing = FALSE, ...)

  \method{csTopTable}{array} (x, ...)

  \method{csTopTable}{NMFfit} (x, ...)

  \method{csTopTable}{NMFfitX} (x, ...)
}
\arguments{
  \item{x}{data object, typically returned by
  \code{\link{ged}}.}

  \item{...}{extra parameters passed to subsequent calls.}

  \item{n}{maximum number of features to extract -- in each
  cell type.}

  \item{ct}{specifies the cell type for which one wants to
  extract the top features.}

  \item{decreasing}{logical that indicates the feature
  ordering, based on their p-values or FDRs.}
}
\description{
  \code{csTopTable} is a generic function that returns
  statistics for each genes, at the cell type levels. It is
  used on deconvolution results from \code{\link{ged}}.
}

