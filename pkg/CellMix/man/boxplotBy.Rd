\name{boxplotBy}
\alias{boxplotBy}
\alias{boxplotBy.default}
\alias{boxplotBy.NMF}
\title{Split Boxplot by Group}
\usage{
  boxplotBy(x, ...)

  \method{boxplotBy}{default} (x, by, scale = FALSE,
    col = rainbow(nlevels(by)), legend = TRUE,
    names = colnames(x), ...)

  \method{boxplotBy}{NMF} (x, by, scale = FALSE, ...)
}
\arguments{
  \item{x}{data matrix, with samples in column.

  NB: in \code{\link{boxplot}} \var{x} has samples in
  rows.}

  \item{by}{grouping variable, coerced to a factor if
  necessary. It is recycled if not of length the number of
  columns in \var{x}}

  \item{col}{spcecifies the colour to use for each group by
  argument \var{by}.}

  \item{legend}{legend specification. If a logical, it
  turns on/off the legend, if a character string, it
  indicates the position of the legend, e.g.
  \code{'topright'} (see \code{\link{legend}} for possible
  values).}

  \item{...}{extra arguments passed to the base
  \code{\link{boxplot}} function.}

  \item{scale}{logical that indicate if each data column
  should be scaled so that it sums up to one.}

  \item{names}{group labels which will be printed under
  each boxplot.  Can be a character vector or an
  \link{expression} (see \link{plotmath}).}
}
\description{
  Plots a boxplot with each data split into sub-boxplots,
  defined according a factor variable.

  The default method plots each \strong{column} using one
  boxplot for each level in a factor.

  The method for \code{\linkS4class{NMF}} objects plots use
  the data from the coefficient matrix.
}
\examples{
x <- rmatrix(3, 20)
boxplotBy(x, gl(2, 10))
}
\author{
  This function was adapted from the following post on
  R-help:

  \url{https://stat.ethz.ch/pipermail/r-help/2003-March/031046.html}
}

