\docType{methods}
\name{markermap}
\alias{basismarkermap}
\alias{markermap}
\alias{markermap,MarkerList,ExpressionSet-method}
\alias{markermap,MarkerList,matrix-method}
\alias{markermap,MarkerList,NMFfitX-method}
\alias{markermap,MarkerList,NMF-method}
\alias{markermap,MatrixData,ANY-method}
\alias{markermap-methods}
\title{Heatmaps Highlighting Markers}
\usage{
  markermap(object, data, ...)

  \S4method{markermap}{MarkerList,matrix}(object, data,
    annCol = NA, annColors = NA, annRow = NA,
    view = "split", subsetRow = NULL, scale = "row",
    color = "YlOrRd:100", Rowv = TRUE, Colv = NA,
    distfun = "correlation", hclustfun = "average", ...)

  basismarkermap(object, data, scale = "r1",
    view = "predict", subsetRow = TRUE, Rowv = NA,
    labRow = NA, ...)
}
\arguments{
  \item{object}{a \code{MarkerList} object}

  \item{data}{reference data object, whose values are used
  in the heatmap.}

  \item{subsetRow}{this argument acts as in
  \code{\link{aheatmap}}, but if \code{subsetRow=TRUE},
  then the heatmap is limited to the markers only.}

  \item{...}{graphical parameters for the text used in
  plot. Parameters passed to \code{\link{grid.text}}, see
  \code{\link{gpar}}.}

  \item{annCol}{specifications of column annotation tracks
  displayed as coloured rows on top of the heatmaps. The
  annotation tracks are drawn from bottom to top. A single
  annotation track can be specified as a single vector;
  multiple tracks are specified as a list, a data frame, or
  an
  \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}
  object, in which case the phenotypic data is used
  (\code{pData(eset)}). Character or integer vectors are
  converted and displayed as factors. Unnamed tracks are
  internally renamed into \code{Xi}, with i being
  incremented for each unamed track, across both column and
  row annotation tracks. For each track, if no
  corresponding colour is specified in argument
  \code{annColors}, a palette or a ramp is automatically
  computed and named after the track's name.}

  \item{annColors}{list for specifying annotation track
  colors manually. It is possible to define the colors for
  only some of the annotations. Check examples for
  details.}

  \item{annRow}{specifications of row annotation tracks
  displayed as coloured columns on the left of the
  heatmaps. The annotation tracks are drawn from left to
  right. The same conversion, renaming and colouring rules
  as for argument \code{annCol} apply.}

  \item{scale}{character indicating how the values should
  scaled in either the row direction or the column
  direction. Note that the scaling is performed after
  row/column clustering, so that it has no effect on the
  row/column ordering. Possible values are: \itemize{ \item
  \code{"row"}: center and standardize each row separately
  to row Z-scores \item \code{"column"}: center and
  standardize each column separately to column Z-scores
  \item \code{"r1"}: scale each row to sum up to one \item
  \code{"c1"}: scale each column to sum up to one \item
  \code{"none"}: no scaling }}

  \item{color}{colour specification for the heatmap.
  Default to palette '-RdYlBu2:100', i.e. reversed palette
  'RdYlBu2' (a slight modification of RColorBrewer's
  palette 'RdYlBu') with 100 colors. Possible values are:
  \itemize{ \item a character/integer vector of length
  greater than 1 that is directly used and assumed to
  contain valid R color specifications. \item a single
  color/integer (between 0 and 8)/other numeric value that
  gives the dominant colors. Numeric values are converted
  into a pallete by \code{rev(sequential_hcl(2, h = x, l =
  c(50, 95)))}. Other values are concatenated with the grey
  colour '#F1F1F1'. \item one of RColorBrewer's palette
  name, or one of 'RdYlBu2', 'rainbow', 'heat', 'topo',
  'terrain', 'cm'. } When the coluor palette is specified
  with a single value, and is negative or preceded a minus
  ('-'), the reversed palette is used. The number of breaks
  can also be specified after a colon (':'). For example,
  the default colour palette is specified as
  '-RdYlBu2:100'.}

  \item{Rowv}{clustering specification(s) for the rows. It
  allows to specify the
  distance/clustering/ordering/display parameters to be
  used for the \emph{rows only}. Possible values are:
  \itemize{ \item \code{TRUE} or \code{NULL} (to be
  consistent with \code{\link{heatmap}}): compute a
  dendrogram from hierarchical clustering using the
  distance and clustering methods \code{distfun} and
  \code{hclustfun}.

  \item \code{NA}: disable any ordering. In this case, and
  if not otherwise specified with argument
  \code{revC=FALSE}, the heatmap shows the input matrix
  with the rows in their original order, with the first row
  on top to the last row at the bottom. Note that this
  differ from the behaviour or \code{\link{heatmap}}, but
  seemed to be a more sensible choice when vizualizing a
  matrix without reordering.

  \item an integer vector of length the number of rows of
  the input matrix (\code{nrow(x)}), that specifies the row
  order. As in the case \code{Rowv=NA}, the ordered matrix
  is shown first row on top, last row at the bottom.

  \item a character vector or a list specifying values to
  use instead of arguments \code{distfun}, \code{hclustfun}
  and \code{reorderfun} when clustering the rows (see the
  respective argument descriptions for a list of accepted
  values). If \code{Rowv} has no names, then the first
  element is used for \code{distfun}, the second (if
  present) is used for \code{hclustfun}, and the third (if
  present) is used for \code{reorderfun}.

  \item a numeric vector of weights, of length the number
  of rows of the input matrix, used to reorder the
  internally computed dendrogram \code{d} by
  \code{reorderfun(d, Rowv)}.

  \item \code{FALSE}: the dendrogram \emph{is} computed
  using methods \code{distfun}, \code{hclustfun}, and
  \code{reorderfun} but is not shown.

  \item a single integer that specifies how many subtrees
  (i.e. clusters) from the computed dendrogram should have
  their root faded out. This can be used to better
  highlight the different clusters.

  \item a single double that specifies how much space is
  used by the computed dendrogram. That is that this value
  is used in place of \code{treeheight}. }}

  \item{Colv}{clustering specification(s) for the columns.
  It accepts the same values as argument \code{Rowv}
  (modulo the expected length for vector specifications),
  and allow specifying the
  distance/clustering/ordering/display parameters to be
  used for the \emph{columns only}. \code{Colv} may also be
  set to \code{"Rowv"}, in which case the dendrogram or
  ordering specifications applied to the rows are also
  applied to the columns. Note that this is allowed only
  for square input matrices, and that the row ordering is
  in this case by default reversed (\code{revC=TRUE}) to
  obtain the diagonal in the standard way (from top-left to
  bottom-right). See argument \code{Rowv} for other
  possible values.}

  \item{distfun}{default distance measure used in
  clustering rows and columns. Possible values are:
  \itemize{ \item \code{"correlation"} and all the
  distances supported by \code{\link{dist}} (e.g.
  \code{"euclidean"}). \item an object of class \code{dist}
  such as returned by \code{\link{dist}} or
  \code{\link{as.dist}}. }}

  \item{hclustfun}{default clustering method used to
  cluster rows and columns. Possible values are: \itemize{
  \item a method name (a character string) supported by
  \code{\link{hclust}} (e.g. \code{'average'}). \item an
  object of class \code{hclust} such as returned by
  \code{\link{hclust}} \item a dendrogram }}

  \item{labRow}{labels for the rows.}

  \item{view}{a character string that specifies how markers
  should be annotated: \describe{
  \item{\dQuote{single}}{markers are shown in a single
  track} \item{\dQuote{split}}{each cell type is shown in a
  separate track} \item{\dQuote{predict}}{one track for
  each column in \code{data} is added, and markers are
  shown on the track associated with the column
  corresponding to its maximum value.} }}
}
\description{
  The function \code{markermap} draws a heatmap of a
  reference expression data (e.g., the expression matrix
  from pure samples or estimated cell-specific signatures),
  where marker are annotated by colored bands on the
  left-hand side of the heatmap.

  \code{basismarkermap} calls \code{markermap} with
  arguments tuned so that: no column reordering is
  performed, the rows are scaled to sum up to one, the
  heatmap only shows marker expression values, and the
  markers are placed according to their most expressing
  estimated signature. This view is meant to help in
  assessing the validity of deconvolution results, when
  known markers are available: ideally marker row
  annotations should be composed of monochrome monoblocks.
}
\details{
  Argument \code{view} controls the way markers are
  annotated. On all views, markers are coloured according
  to their type, defined by the element of \code{object} in
  which they appear. Each couloured tick/cell corresponds
  to a different marker position.

  The following views are available: \describe{
  \item{single}{ a single row annotation is added, showing
  the position of each marker.} \item{split}{ one track per
  marker type in added, showing the position of each marker
  in its respective cell type.} \item{predict}{ one track
  per column in \code{data} is added, showing the position
  of each marker in the most expressing column. When
  \code{data} is a basis matrix obtained from
  deconvolution, this view is useful to check how known
  cell type markers (the coulours) map on estimated
  signatures (the annotation columns).} }
}
\section{Methods}{
  \describe{

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "matrix")}: Workhorse method is for
  \code{markermap}. }

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "ExpressionSet")}: The method \code{markermap} for
  \code{ExpressionSet} objects calls the main
  \code{markermap} method on the expression matrix
  \code{exprs(object)}. }

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "NMF")}: The method \code{markermap} for
  \code{NMF} objects calls the main \code{markermap} method
  on the basis matrix \code{basis(object)}. }

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "NMFfitX")}: The method \code{markermap} for
  \code{NMFfitX} objects calls the main \code{markermap}
  method on the best fit \code{fit(object)}. }

  \item{markermap}{\code{signature(object = "MatrixData",
  data = "ANY")}: This method extracts and plots a list of
  markers from a set of basis signature matrix (in
  \code{object}), where each feature is associated with the
  most-expressing signature. }

  }
}
\examples{
x <- rmix(3, 100, 20)
m <- getMarkers(x)
markermap(m, basis(x))
markermap(m, x, view='single')
basismarkermap(m, rnmf(3, x))

# after real deconvolution
res <- ged(x, coef(x), 'csSAM')
basismarkermap(m, res)
markermap(m, res, view='split')
}
\keyword{methods}

