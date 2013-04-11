\name{applyBy}
\alias{applyBy}
\alias{colApplyBy}
\alias{colMaxsBy}
\alias{colMeansBy}
\alias{colMediansBy}
\alias{colSumsBy}
\alias{rowApplyBy}
\alias{rowMaxsBy}
\alias{rowMeansBy}
\alias{rowMediansBy}
\alias{rowSumsBy}
\title{Group Apply}
\usage{
  applyBy(x, BY, MARGIN, FUN, W = NULL, ..., DROP = FALSE)

  rowApplyBy(x, BY, FUN, ...)

  colApplyBy(x, BY, FUN, ...)

  colSumsBy(x, BY, ...)

  rowSumsBy(x, BY, ...)

  rowMeansBy(x, BY, ...)

  colMeansBy(x, BY, ...)

  rowMediansBy(x, BY, ...)

  colMediansBy(x, BY, ...)

  rowMaxsBy(x, BY, ...)

  colMaxsBy(x, BY, ...)
}
\arguments{
  \item{x}{matrix-like object on which \code{\link{apply}}
  can be called.}

  \item{BY}{factor or object coerced to a factor, that
  defines the groups within which the function \code{FUN}
  is applied.}

  \item{MARGIN}{margin along which the function \code{FUN}
  is applied: 1L for rows, 2L for columns.}

  \item{FUN}{function to apply to each sub-matrix that
  contains the rows/columns defined by each level of
  argument \code{BY}. It must be a function that takes a
  matrix as its first argument and returns a vector of
  length the dimension of margin \code{MARGIN} of
  \code{x}.}

  \item{...}{extra parameters passed to \code{FUN}.}

  \item{DROP}{logical that indicates if absent levels
  should be removed from the result matrix, or appear as
  0-filled rows/columns.}

  \item{W}{An optional \code{\link[base]{numeric}} NxM
  \code{\link[base]{matrix}} of weights.}
}
\value{
  The result is a matrix whose margin's dimension
  \code{MARGIN} is equal the same margin's dimension in
  \code{x}, and the other to the number of levels in
  \code{BY}.
}
\description{
  Applies a function to rows or columns of a matrix-like
  object, separately within groups defined by a factor.

  \code{rowApplyBy} applies a function to rows of
  sub-matrices whose columns are defined by a factor.

  \code{rowApplyBy} applies a function to columns of
  sub-matrices whose rows are defined by a factor.

  \code{col<STAT>By} computes for each column a given
  statistic within separate groups of rows, which are
  defined by a factor.

  \code{row<STAT>By} computes for each row a given
  statistic within separate groups of columns, which are
  defined by a factor.
}
\details{
  \code{applyBy} is the work horse function that is called
  by other more user-friendly functions.

  \code{applyBy} is a wrapper around
  \code{\link[matrixStats]{colAvgsPerRowSet}}, which make
  the computation really fast, but requires somehow
  cumbersome matrix specifications for the groups of
  columns or rows. The wrapper builds the arguments for the
  particular case where the groups are defined by a factor.
}
\examples{
# random data matrix
x <- rmatrix(12, 6)

# by groups of columns
fc <- gl(2, 3)
b <- applyBy(x, fc, 1L, rowSums)
b
# or
balt <- rowApplyBy(x, fc, rowSums)
identical(b, balt)

# by groups of rows
fr <- gl(3, 4)
b <- applyBy(x, fr, 2L, colSums)
# or
balt <- colApplyBy(x, fr, colSums)
identical(b, balt)
}

