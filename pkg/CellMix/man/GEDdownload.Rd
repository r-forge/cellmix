\name{GEDdownload}
\alias{GEDdownload}
\title{Downloading Gene Expression Deconvolution Datasets}
\usage{
  GEDdownload(x, ..., destdir = GEDtmp(),
    datasource = NULL, annotation = NULL, verbose = TRUE)
}
\arguments{
  \item{x}{access key of a registered dataset}

  \item{...}{extra arguments passed to
  \code{\link{download.file}} and/or \code{\link{getGSE}}.}

  \item{destdir}{destination directory for all the
  downloaded files}

  \item{datasource}{Data source from where to fetch the
  dataset, e.g., an online dataset repository lke
  \emph{GEO} or \emph{ArrayExpress}. If missing, it is
  determined automatically from the access key, but is
  honoured if provided. Currently only values \code{'GEO'}
  and \code{'ArrayExpress'} (or its alias \code{'AE'}) are
  supported.}

  \item{annotation}{annotation package (a string) to attach
  to the created ExpressionSet object.}

  \item{verbose}{logical or numeric that indicate the
  verbosity level}
}
\value{
  the value returned by \code{\link{getGSE}}.
}
\description{
  \code{GEDdownload} downloads gene expression datasets
  that for datasets that have an entry in the \pkg{CellMix}
  package's internal registry.
}

