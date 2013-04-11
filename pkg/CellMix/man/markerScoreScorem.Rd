\name{markerScoreScorem}
\alias{markerScoreScorem}
\alias{selectMarkers.markerScore_scorem}
\title{Marker scoring method: SCOREM}
\usage{
  markerScoreScorem(object, data, minsd = 0.1,
    alpha = 0.01, ..., verbose = FALSE)

  \method{selectMarkers}{markerScore_scorem} (x, data,
    mscore = "pearson", statistic = c("size", "score"),
    ..., .object = NULL)
}
\arguments{
  \item{minsd}{minimum standard deviation of expression
  marker profiles}

  \item{alpha}{numerical value between 0 and 1, indicating
  the level of significance required for correlations to be
  considered significant. It is a p-value threshold: the
  lower the value the more markers must be correlated to be
  called in the same sub-group.}

  \item{...}{extra arguments not used in
  \code{markerScoreScorem}, but passed to
  \code{\link{selectMarkers.MarkerList}} in
  \code{selectMarkers}.}

  \item{verbose}{logical that toggles log messages (passed
  to \code{\link{subset,MarkerList-method}}), which can be
  useful to track any ID conversion.}

  \item{mscore}{correlation method used to compute a score
  for each marker of the selected groups.}

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

  \item{x}{data object based on which the markers are
  selected, as computed returned by
  \code{extractMarkers(..., format='raw')}. The type of
  \code{x} depends on the scoring method used to compute
  it.}

  \item{statistic}{name or column index of the
  statistic/score to use. The default is to use the second
  column of the matrix returned by the scoring method, but
  each method can define its own default.}

  \item{.object}{argument used internally by
  \code{\link{extractMarkers}} to pass the original object
  on which marker scores where computed.}
}
\description{
  This method applies an approach called \emph{statistical
  consolidation of redundant expression measures} (SCOREM)
  from \cite{Schneider et al. (2011)} whitin each cell
  type, in order to select groups of markers that have
  consistent expression in mixed samples.

  If multiple groups pass the criterion, then the group
  with the lowest median gene expression gene-wise standard
  deviation is selected. This is to select the group of
  markers whose expression values are the more similar.
}
\details{
  The SCOREM approach allows sets of marker genes to be
  optimised with respect to the dataset one intends to use
  them with. It uses a sub-graph detection algorithm, on a
  graph defined from the rank correlation matrix of marker
  expression profiles.

  \strong{IMPORTANT: this feature is still experimental and
  in development.}

  Before returning the \code{MarkerList} object, each
  selected marker is assigned a score defined as its
  minimum correlation with other markers in its group.
}
\examples{
# random data matrix and marker list
x <- rmix(3, 100, 20)
m <- getMarkers(x)

# compute SCOREM groups: should keep all markers together
# or remove weak markers
s <- markerScoreScorem(m, x)
str(s)
# add some non marker genes
m2 <- m
m2[[1]] <- as.integer(c(m2[[1]], c(40, 50, 60)))
s <- markerScoreScorem(m2, x)
str(s[[1]])

# with group selection and reordering through the main interface function
extractMarkers(m, x, method='scorem')
extractMarkers(m2, x, method='scorem')
}
\references{
  Schneider S, Smith T and Hansen U (2011). "SCOREM:
  statistical consolidation of redundant expression
  measures." _Nucleic acids research_, pp. 1-12. ISSN
  1362-4962, <URL: http://dx.doi.org/10.1093/nar/gkr1270>,
  <URL: http://www.ncbi.nlm.nih.gov/pubmed/22210887>.
}
\seealso{
  Other markerScore: \code{\link{extractMarkers}},
  \code{\link{markerScoreAbbas}},
  \code{\link{markerScoreHSD}},
  \code{\link{markerScoreMaxcol}},
  \code{\link{markerScoreMethod}},
  \code{\link{scoreMarkers}},
  \code{\link{selectMarkers.markerScore_HSD}}
}

