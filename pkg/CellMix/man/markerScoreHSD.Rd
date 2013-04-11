\name{markerScoreHSD}
\alias{markerScoreHSD}
\alias{selectMarkers.markerScore_HSD}
\title{Marker Scoring Method: Tukey Honest Significant Difference}
\usage{
  markerScoreHSD(object, data, log = !is_logscale(object),
    lbase = 2, verbose = FALSE)

  \method{selectMarkers}{markerScore_HSD} (x, data,
    statistic = max, ...)
}
\arguments{
  \item{verbose}{verbosity level, usually \code{TRUE} or
  \code{FALSE}.}

  \item{statistic}{method used to aggregate p-values of the
  pairwise comparisons with each of the other cell types
  into a single numeric score.}

  \item{object}{A numeric \code{matrix}, an object of class
  \code{\link[Biobase]{ExpressionSet}}, or a
  \code{\linkS4class{MarkerList}} object.}

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

  \item{log}{a logical specifying if the should the data be
  log-transformed before computing the p-values. Default is
  to log-transform the data if is not already in log-scale,
  which is determined by the function
  \code{\link{is_logscale}}.}

  \item{lbase}{log base to used if the data is
  log-transformed.}

  \item{x}{data object based on which the markers are
  selected, as computed returned by
  \code{extractMarkers(..., format='raw')}. The type of
  \code{x} depends on the scoring method used to compute
  it.}

  \item{...}{other parameters passed to the scoring method
  in \code{extractMarkers} \strong{and}
  \code{\link{selectMarkers}}, or that define default
  arguments when defining a scoring method with
  \code{markerScoreMethod}.}
}
\description{
  This scoring method \code{markerScoreHSD} performs
  pairwise comparison between groups of pure samples, and
  scores each comparison using Tukey's Honest Significant
  Difference p-values (see \code{\link{TukeyHSD}}).

  The method \code{selectMarkersmarkerScore_HSD} selects,
  within each cell type separately, the markers with the
  lower aggregated p-value for Tukey's HSD. The default
  aggregation method is to compute the maximum HSD p-value.
}
\details{
  The scores are returned in a matrix, with features in
  rows and cell types in column, which contains the HSD
  p-values corresponding to the comparisons between the
  most expressing cell type and other cell types. Each row
  contains an NA value that identifies the column
  corresponding to the associated feature's most expressing
  cell type.

  Features whose expression is not consistently higher in
  one cell type than in any other cell type are discarded.
}
\examples{
# generate data from pure cell type samples
x <- rpure(3)
x
aheatmap(x, annCol=TRUE)

# extract markers
ml <- extractMarkers(x, x$CellType, method='HSD')
# check score/p-value distribution
hist(ml)
# plot most significant ones
profplot(ml[ml < 0.0001], x, split=TRUE)
}
\seealso{
  \code{\link{TukeyHSD}}

  Other markerScore: \code{\link{extractMarkers}},
  \code{\link{markerScoreAbbas}},
  \code{\link{markerScoreMaxcol}},
  \code{\link{markerScoreMethod}},
  \code{\link{markerScoreScorem}},
  \code{\link{scoreMarkers}},
  \code{\link{selectMarkers.markerScore_scorem}}
}

