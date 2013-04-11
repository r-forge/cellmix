\docType{methods}
\name{.atrack,MarkerList-method}
\alias{.atrack,MarkerList-method}
\title{Heatmap Annotation Track for MarkerList Objects}
\usage{
  \S4method{.atrack}{MarkerList}(object, data = NULL, ...,
    view = c("single", "split", "predict"))
}
\arguments{
  \item{object}{a MarkerList object}

  \item{data}{the reference data to annotate}

  \item{view}{a character string that specifies how markers
  should be annotated: \describe{
  \item{\dQuote{single}}{markers are shown in a single
  track} \item{\dQuote{split}}{each cell type is shown in a
  separate track} \item{\dQuote{predict}}{one track for
  each column in \code{data} is added, and markers are
  shown on the track associated with the column
  corresponding to its maximum value.} }}

  \item{...}{extra parameters passed to .atrack}
}
\description{
  This method allows to add annotation tracks in heatmaps
  produced by \code{\link{aheatmap}}, to highlight the
  position of markers, when plotting either the global
  expression values or cell type-specific signatures,
  whether measured or estimated.
}
\examples{
# load IRIS markers and the Abbas signature matrix
m <- MarkerList('IRIS')
data(Abbas)
aheatmap(Abbas, annRow=m)
}
\keyword{methods}

