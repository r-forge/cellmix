\name{csplot}
\alias{csplot}
\alias{csplot.NMFfit}
\alias{csplot.NMFfitX}
\title{Plots Cell-Specific FDR Estimates}
\usage{
  csplot(x, ..., filename = NULL)

  \method{csplot}{NMFfit} (x, ...)

  \method{csplot}{NMFfitX} (x, ...)
}
\arguments{
  \item{x}{data object, typically returned by
  \code{\link{ged}}.}

  \item{...}{extra parameters passed to subsequent calls.}

  \item{filename}{file where to save the plot. Can be of
  any type amongst pdf, svg, png, bmp, jpg, tiff.}
}
\description{
  Plots Cell-Specific FDR Estimates
}

