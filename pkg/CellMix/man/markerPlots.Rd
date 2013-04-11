\name{beeswarm.MarkerList}
\alias{beeswarm.MarkerList}
\alias{hist.MarkerList}
\alias{profplot.MarkerList}
\alias{screeplot.MarkerList}
\alias{stripchart.MarkerList}
\title{Plotting Markers}
\usage{
  \method{beeswarm}{MarkerList} (x, data, types = NULL,
    tcol = "red", otcol = "#00000090", las = 2,
    fold = FALSE, method = "square", corral = "wrap",
    pch = 19, pwcol = NULL,
    ylab = if (fold) "Fold changes" else "Expression values",
    ...)

  \method{stripchart}{MarkerList} (x, data, types = NULL,
    tcol = "red", otcol = "#00000090", las = 2,
    fold = FALSE, method = "square", corral = "wrap",
    pch = 19, pwcol = NULL,
    ylab = if (fold) "Fold changes" else "Expression values",
    ...)

  \method{profplot}{MarkerList} (x, y, groups = NULL,
    scale = FALSE, col = NULL, legend = FALSE,
    split = TRUE, ylab = NA, labels = NULL, ...,
    restore.gpar = TRUE)

  \method{screeplot}{MarkerList} (x, data = NULL,
    breakdown = FALSE, range = NULL, xlab = NULL,
    ylab = NULL, main = NULL, ..., plot = TRUE)

  \method{hist}{MarkerList} (x, range = NULL,
    split = FALSE, xlab = "Marker scores",
    main = "Histogram of marker scores", ...,
    restore.gpar = TRUE)
}
\arguments{
  \item{x}{an object of class \code{\link{MarkerList}}}

  \item{data}{a numeric \code{matrix} or an
  \code{ExpressionSet} object containing expression values
  from samples.}

  \item{types}{a \code{factor} giving the cell-type of each
  sample}

  \item{tcol}{type colour, used for the marker values in
  samples of their corresponding types.}

  \item{otcol}{out-of-type colour, used for maker values in
  samples that are not from their corresponding type.}

  \item{las}{specifies the orientation of the y-labels: 1
  for horizontal, 2 for vertical. See \code{\link{par}}.}

  \item{fold}{logical that indicates if the plot should
  show the markers' plain expression values or their fold
  change against the maximum expression value in samples
  from other types.}

  \item{...}{extra parameters passed to subsequent calls to
  \code{\link{beeswarm}}, \code{\link{matplot}},
  \code{\link{barplot}}, \code{\link{hist}} or
  \code{\link{plot}}, depending on the main calling
  function.}

  \item{y}{reference data matrix, that contains the marker
  expression values that are plotted.}

  \item{groups}{factor coercible vector that defines groups
  of samples, which are plotted using different point
  types.}

  \item{scale}{when \code{split=FALSE}, this is a logical
  that indicates if the marker expression profiles should
  be scaled into relative contribution (i.e., sum up to
  one), as in \code{\link{profplot}}.

  When \var{split} is not \code{FALSE}, then
  \code{scale=TRUE} indicates that the profiles should be
  normalised with the mean expression profile of all
  markers. If a single numeric, then it indicates the index
  of a specific marker that is used to normalise the
  expression profiles of other markers. Separate indexes
  for each cell type can be passed as a numeric vector or a
  list.

  This is useful to highlight groups of markers that should
  have equivalent cell-specific expression levels, given
  the -- mixed -- expression data in \code{y}. These groups
  are likely to provide better proportion/cell-specific
  estimates.

  Note: values are recycled if necessary.}

  \item{col}{colour specification for each cell types if
  \code{split=FALSE}, or for each marker in each cell type
  otherwise.}

  \item{split}{logical that indicates if one should draw a
  plot for each type separately.

  One can specify how to split the plots, by providing a
  2-long vector that is used as \code{par(mfrow=split)}.}

  \item{ylab}{label for the y-axis. If missing, a default
  label is generated, which, in particular, tells if the
  plot is scaled or not.}

  \item{restore.gpar}{logical that indicates if all
  graphical parameters should be restored on exit. It could
  be used to add things to the plots, especially when
  \code{split=TRUE}.}

  \item{breakdown}{logical that indicates if breakdowns of
  the number of markers over the given range of values
  should be shown. It may also be a single numberic that
  indicates the number of breakdowns to show.}

  \item{range}{specifies the range of values to plot. In
  \code{hist}, or in \code{screeplot} when
  \code{data=NULL}, it indicates a range over which the
  scores are used to computing the histogram. If
  \code{NULL}, then all values are used (Note: this might
  take a bit of time for long marker lists if
  \code{breakdown=TRUE}). Otherwise it must be a 2-long
  numeric vector giving the interval of values to use.

  If \code{data} is not \code{NULL}, it indicates which}

  \item{plot}{logical that indicates if the function should
  effectively produce a plot or only return the result.

  Otherwise, if \var{data} is not \code{NULL}, then it must
  be either a range of values}

  \item{method}{ Method for arranging points (see Details).
  }

  \item{corral}{ Method to adjust points that would be
  placed outside their own group region (see Details). }

  \item{pch}{ Plotting characters and colors, specified by
  group (see Details). }

  \item{pwcol}{ \dQuote{Point-wise} plotting characters and
  colors, specified for each data point (see Details). }

  \item{labels}{ Labels for each group.  By default, these
  are inferred from the data. }

  \item{xlab}{ Axis labels. }

  \item{legend}{a logical that specifies whether drawing
  the legend or not, or coordinates specifications passed
  to argument \code{x} of \code{\link{legend}}, that
  specifies the position of the legend.}

  \item{main}{graphics parameters.}
}
\value{
  When \code{data=NULL} and \code{breakdown} is not
  \code{FALSE} this function returns a list with the result
  from the breakdown barplot and the breakdown (in element
  \sQuote{breakdown}).

  When \code{data} is not \code{NULL}, a
  \code{\linkS4class{MarkerList}} object made of the set of
  markers that achieve the least condition number is
  returned.
}
\description{
  The method \code{beeswarm} for \code{MarkerList} objects
  draws, separately for each cell-type, a stripchart of the
  markers' expression values in a reference data,
  highlighting the expression of their respective
  associated -- generally pure -- samples.

  The method \code{stripchart.MarkerList} is an alias for
  \code{beeswarm.MarkerList}.

  \code{profplot} is vectorised over arguments \var{scale}
  and \var{legend}.

  \code{screeplot} plots the following diagnostic plots:
  \itemize{ \item for marker lists with associated scores,
  it plots the total number of markers whose values are
  over/under a range of threshold, as well as, optionally
  their distribtuion amongst the different cell types;
  \item a screeplot of the condition number of a
  deconvolution matrix for a varying (increasing) number of
  markers. }

  \code{hist} plots the histogram of numeric scores
  associated with each marker in a marker list, if any.
}

