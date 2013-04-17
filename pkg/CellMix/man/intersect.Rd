\docType{methods}
\name{intersect}
\alias{intersect}
\alias{intersect,MatrixData,character-method}
\alias{intersect,MatrixData,logical-method}
\alias{intersect,MatrixData,MarkerList-method}
\alias{intersect,MatrixData,MatrixData-method}
\title{Enhanced Subsetting for Matrix-like Data}
\usage{
  \S4method{intersect}{MatrixData,logical}(x, y)

  \S4method{intersect}{MatrixData,character}(x, y)

  \S4method{intersect}{MatrixData,MatrixData}(x, y)

  \S4method{intersect}{MatrixData,MarkerList}(x, y)
}
\arguments{
  \item{x}{vectors (of the same mode) containing a sequence
  of items (conceptually) with no duplicated values.}

  \item{y}{vectors (of the same mode) containing a sequence
  of items (conceptually) with no duplicated values.}
}
\description{
  These methods subset Matrix-like data only keeping rows
  whose names are in common with some other data.
}
\details{
  The generic for these methods is taken from the
  \pkg{GSEABase} package because it is entirely imported by
  \pkg{CellMix}. However, a more natural definition would
  be to be able to use the -- currently conflicting --
  version defined in the recommended \pkg{BiocGenerics}
  package.
}
\section{Methods}{
  \describe{

  \item{intersect}{\code{signature(x = "MatrixData", y =
  "logical")}: This method is equivalent to
  \code{\link[base]{subset}}. }

  \item{intersect}{\code{signature(x = "MatrixData", y =
  "character")}: Subset a matrix-like object by only
  keeping the rows whose names are in a given reference
  character vector. }

  \item{intersect}{\code{signature(x = "MatrixData", y =
  "MatrixData")}: Subset a matrix-like object by only
  keeping the rows whose names are in common with another
  matrix-like data.

  This is a shortcut for \code{intersect(x,
  featureNames(y), ...)} }

  \item{intersect}{\code{signature(x = "MatrixData", y =
  "MarkerList")}: Subset a matrix-like object by only
  keeping the features in common with a \code{MarkerList}
  object. Note that no ID conversion is attempted, meaning
  that marker IDs must match exactly feature names used in
  the matrix data.

  This is a shortcut for \code{intersect(x,
  featureNames(y), ...)}. }

  }
}
\keyword{methods}

