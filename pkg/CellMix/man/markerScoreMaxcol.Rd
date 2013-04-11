\name{markerScoreMaxcol}
\alias{markerScoreMaxcol}
\title{Marker scoring method: Max Expression}
\usage{
  markerScoreMaxcol(object)
}
\arguments{
  \item{object}{A numeric \code{matrix}, an object of class
  \code{\link[Biobase]{ExpressionSet}}, or a
  \code{\linkS4class{MarkerList}} object.}
}
\description{
  This method scores genes according with their highest
  expression in a set of pure samples.
}
\details{
  Genes are assigned to their respective maximum expressing
  column/value and returned as a \code{MarkerList} object.
}
\examples{
# random data matrix with rownames (needed)
x <- rmatrix(20, 3)
rownames(x) <- paste("g", 1:nrow(x), sep='_')

# compute score
markerScoreMaxcol(x)
# same thing with auto-ordering through the main interface function
extractMarkers(x, method='maxcol')
}
\seealso{
  Other markerScore: \code{\link{extractMarkers}},
  \code{\link{markerScoreAbbas}},
  \code{\link{markerScoreHSD}},
  \code{\link{markerScoreMethod}},
  \code{\link{markerScoreScorem}},
  \code{\link{scoreMarkers}},
  \code{\link{selectMarkers.markerScore_HSD}},
  \code{\link{selectMarkers.markerScore_scorem}}
}

