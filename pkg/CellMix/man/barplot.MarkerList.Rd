\name{barplot.MarkerList}
\alias{barplot.MarkerList}
\title{Plotting MarkerList Objects}
\usage{
  \method{barplot}{MarkerList} (height, data = NULL,
    col = NULL, byType = TRUE, ylab = NULL,
    border = missing(data), beside = TRUE, scale = !beside,
    legend = TRUE,
    reorder = if (is.null(data)) -1L else FALSE, las = 2,
    ...)
}
\arguments{
  \item{height}{a \code{\linkS4class{MarkerList}} object.}

  \item{data}{dataset (a matrix or
  \code{\link{ExpressionSet}} object) that contains
  expression values of the markers in some experiment. If
  not missing, values corresponding to markers are plotted
  for each column in the data. Markers that do not have
  corresponding values in \code{data} are discarded. Note
  that no identifier conversion is performed, meaning that
  the marker list and the dataset must use the type of gene
  identifiers. Use the function \code{\link{convertIDs}},
  if a conversion is needed.}

  \item{col}{colour specification for the bars.}

  \item{byType}{logical that indicates if the colour
  specification in \code{col} should be interpreted as
  colours for each cell type. Is is used only when
  \code{data} is not missing. If \code{TRUE} then
  \code{col} specifies the colours to be used for all bars
  in each cell type. If \code{FALSE}}

  \item{ylab}{label for the y-axis. Default labels are
  provided if missing.}

  \item{scale}{logical that indicates if the data should be
  scaled. This argument is used only when
  \code{beside=FALSE}. If \code{TRUE} the gene expression
  profiles (rows of \code{data}) are scaled to sum up to
  100, then the sample expression profiles (the columns of
  \code{data}) are scaled to sum up to 100. The scaled
  values can be interpreted as how much a marker gene is
  expressed in a sample relatively to other samples and
  genes.}

  \item{legend}{logical that indicates if a legend box
  should be added.}

  \item{reorder}{when \code{data} is missing, this argument
  indicates if/how the cell types should be reordered:
  \code{FALSE} do not reorder the cell types, a negative
  number orders by decreasing number of markers, and
  anything other than \code{FALSE} orders by increasing
  order of number of markers.

  When \code{data} is not missing, then it indicates if the
  markers should be ordered in the same order as they
  appear in the data.

  is used when \code{data} is not missing and indicates if
  the markers should be ordered as they appear in the
  dataset, as opposed to keep their original order, i.e.
  grouped by cell type.}

  \item{las}{specifies the orientation of the x-labels: 1
  for horizontal, 2 for vertical.}

  \item{...}{extra arguments passed to
  \code{\link[graphics]{barplot}}.}

  \item{border}{the color to be used for the border of the
  bars.  Use \code{border = NA} to omit borders.  If there
  are shading lines, \code{border = TRUE} means use the
  same colour for the border as for the shading lines.}

  \item{beside}{a logical value.  If \code{FALSE}, the
  columns of \code{height} are portrayed as stacked bars,
  and if \code{TRUE} the columns are portrayed as
  juxtaposed bars.}
}
\description{
  The function \code{barplot} plots a barplot of the number
  of markers of each cell type present in a
  \code{\linkS4class{MarkerList}} object, or their
  respective expression values in a given dataset.
}

