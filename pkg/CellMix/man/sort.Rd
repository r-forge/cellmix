\docType{methods}
\name{reorder,MarkerList-method}
\alias{reorder,MarkerList-method}
\alias{sort.MarkerList}
\title{Reordering Marker Lists}
\usage{
  \S4method{reorder}{MarkerList}(x, data, by = NULL,
    fun = max, by.fun = rowMins,
    decreasing = !missing(data))

  \method{sort}{MarkerList} (x, decreasing = FALSE, ...)
}
\arguments{
  \item{x}{a \code{MarkerList} object}

  \item{data}{reference data as a numeric vector, a matrix,
  or an \code{ExpressionSet} object. These data must have
  names or rownames against which markers can be match.}

  \item{by}{optional \code{factor} only used when
  \code{data} is a matrix-like object, to}

  \item{fun}{function to perform the second step in the
  computation of scores (cf. argument \code{by}).}

  \item{by.fun}{function to perform the first step in the
  computation of scores (cf. argument \code{by}).}

  \item{decreasing}{logical that indicates how the markers
  should be ordered according to their score.}

  \item{...}{ optional: extra arguments supplied to
  \code{FUN}}
}
\value{
  a \code{MarkerList} object that also contain the computed
  scores.
}
\description{
  The method \code{reorder} for \code{MarkerList} objects
  allows to reorder each set of markers according to some
  auxiliary reference data, e.g, their gene expression in
  different cell types.

  \code{sort} sorts each element of a marker list.
}
\details{
  When \code{data} is a matrix-like object, a score is
  computed for each marker in two steps: the first step
  computes scores for each marker within each group of
  columns defined by the levels of \code{by} (function
  \code{rowMins}); the second step aggregates these scores
  into a single value for each marker (function
  \code{fun}).

  The markers are then ordered wihtin their respective set,
  according to the score and argument \code{decreasing}.
}
\seealso{
  \code{\link{reverse,MarkerList-method}}
}
\keyword{methods}

