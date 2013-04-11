\docType{methods}
\name{extractMarkers}
\alias{extractMarkers}
\alias{extractMarkers,ANY-method}
\alias{extractMarkers,ExpressionSet-method}
\alias{extractMarkers-methods}
\alias{markerScoreMethod}
\alias{scoreMarkers}
\title{Extract Markers from Pure Samples.}
\usage{
  extractMarkers(object, ...)

  \S4method{extractMarkers}{ANY}(object, data = NULL,
    method = "Abbas", threshold = NA, decreasing = FALSE,
    ..., format = c("list", "raw"))

  scoreMarkers(...)

  markerScoreMethod(method, object, ...)
}
\arguments{
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

  \item{method}{accession key of the method either used to
  compute marker scores in \code{extractMarkers}, or to
  define/retrieve scoring methods in
  \code{markerScoreMethod}.

  All available methods can be retrieved via
  \code{markerScoreMethod()}.

  A custom scoring method can also be passed as a
  function.}

  \item{...}{other parameters passed to the scoring method
  in \code{extractMarkers} \strong{and}
  \code{\link{selectMarkers}}, or that define default
  arguments when defining a scoring method with
  \code{markerScoreMethod}.}

  \item{format}{specifies the output format: \sQuote{list}
  returns a \code{\linkS4class{MarkerList}} object and
  \sQuote{raw} directly returns the result of the scoring
  method -- which can be useful for lengthy computations.}

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
}
\description{
  When expression data from pure samples are available, it
  is possible to select genes for their ability to
  discriminate between the different cell-types. The
  function \code{extractMarkers} provides an interface to
  compute scores or statistics using different methods, and
  select the ones that pass a given threshold.

  For long-runnning scoring methods, the marker selection
  process could be performed manually, by separating the
  scoring and selections steps, e.g.: \samp{ sc <-
  scoreMarkers(eset, data=groups, method='HSD')
  selectMarkers(sc, threshold=10^-3) }

  \code{markerScoreMethod} defines/retrieves marker scoring
  methods.
}
\details{
  To see all available methods use
  \code{markerScoreMethod()}.
}
\section{Methods}{
  \describe{

  \item{extractMarkers}{\code{signature(object =
  "ExpressionSet")}: The method \code{extractMarkers} for
  \code{ExpressionSet} objects calls \code{extractMarkers}
  on the expression matrix \code{exprs(object)}. }

  }
}
\examples{
# see available scoring methods
markerScoreMethod()
}
\seealso{
  Other markerScore: \code{\link{markerScoreAbbas}},
  \code{\link{markerScoreHSD}},
  \code{\link{markerScoreMaxcol}},
  \code{\link{markerScoreScorem}},
  \code{\link{selectMarkers.markerScore_HSD}},
  \code{\link{selectMarkers.markerScore_scorem}}
}
\keyword{methods}

