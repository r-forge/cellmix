\name{propplot}
\alias{propplot}
\title{Cell Type Proportion Plot}
\usage{
  propplot(x, col = NULL, legend.text = TRUE,
    scale = FALSE, legwidth = 10, ...)
}
\arguments{
  \item{x}{data. Can be a matrix with cell types in rows
  and samples in columns, or an \code{\linkS4class{NMF}}
  object that contains proportions as its mixture
  coefficient matrix, e.g. as computed by
  \code{\link{ged}}.}

  \item{col}{color specification}

  \item{legend.text}{logical that indicates if a legend
  should be added to the plot, or the text of the legend as
  in \code{\link{barplot}}.}

  \item{scale}{logical that indicates if the data should be
  scaled so that columns sum up to one, like proportions.}

  \item{legwidth}{number of margin lines left to display
  the legend in the plot margin.}

  \item{...}{extra arguments passed to
  \code{\link{barplot}}.}
}
\description{
  Draws a vertical barplot of cell-types proportions.
}

