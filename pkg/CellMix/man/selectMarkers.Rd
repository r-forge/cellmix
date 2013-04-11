\name{selectMarkers}
\alias{selectMarkers}
\alias{selectMarkers.MarkerList}
\alias{selectMarkers.markerScore}
\title{Select Markers Based on Scores}
\usage{
  selectMarkers(x, ..., .object = NULL)

  \method{selectMarkers}{MarkerList} (x, threshold = NA,
    decreasing = FALSE, ..., .object)

  \method{selectMarkers}{markerScore} (x, statistic = 2L,
    data = attr(x, "types"), ..., .object)
}
\arguments{
  \item{x}{data object based on which the markers are
  selected, as computed returned by
  \code{extractMarkers(..., format='raw')}. The type of
  \code{x} depends on the scoring method used to compute
  it.}

  \item{.object}{argument used internally by
  \code{\link{extractMarkers}} to pass the original object
  on which marker scores where computed.}

  \item{threshold}{threshold that is applied to filter
  markers based on their respective statistic/score: Genes
  with \code{statistic >= threshold} are selected if
  \code{decreasing=TRUE}. Otherwise, if
  \code{decreasing=FALSE} or \code{NA}, the selected genes
  are those with \code{statistic <= threshold}.

  Filtering is disabled if \code{threshold=NA}.}

  \item{decreasing}{logical that indicates how the
  statistic/score should be ordered: \code{TRUE} orders by
  decreasing value (i.e. the greater the score the better),
  while \code{FALSE} orders by increasing value (i.e. the
  lower the score the better). The value of this argument
  also affects the way the filtering is performed (see
  argument \code{threshold}).

  Ordering is disabled if \code{decreasing=NA}.}

  \item{statistic}{name or column index of the
  statistic/score to use. The default is to use the second
  column of the matrix returned by the scoring method, but
  each method can define its own default.}

  \item{...}{other parameters passed to the scoring method
  in \code{extractMarkers} \strong{and}
  \code{\link{selectMarkers}}, or that define default
  arguments when defining a scoring method with
  \code{markerScoreMethod}.}

  \item{data}{main extra data used by the scoring and
  selection methods.

  If \var{object} is a matrix-like object, \var{data} is
  generally a factor or grouping variable, as a vector,
  that defines the cell-type for each -- pure -- sample. If
  a vector, this argument is converted into a factor with
  the levels in their order of appearance in \code{data},
  by \code{factor(data, levels=unique(data))}. This is to
  obtain levels in an order that is consistent with the
  samples' order.

  If \var{object} is a \code{MarkerList} object, then
  \var{data} is generally a matrix-like object that
  contains expression data.}
}
\value{
  \code{selectMarkers} returns an object of class
  \code{MarkerList}.
}
\description{
  The function \code{selectMarkers} filters raw results
  returned by scoring methods and wrap them into a
  \code{MarkerList} object. It is called by
  \code{\link{extractMarkers}} when \var{format} is not
  \code{'raw'}.
}
\details{
  The default method \code{selectMarkers.markerScore}
  selects markers based on scores stored as a matrix-like
  object, whose first column contains a factor or a vector
  that assigns each feature to a group (cell type) and
  column with index \code{statistic} (by default the second
  column) contains the feature scores.

  Scoring methods that do not define a dedicated
  \code{selectMarkers} must return a \code{MarkerList}
  object or a matrix that complies with the description
  above. It is recommended that, if such methods compute a
  range of statistics, they store the requested selection
  statistic in the second column.
}

