\name{setGEDMethod}
\alias{removeGEDMethod}
\alias{setGEDMethod}
\title{Register CellMix Deconvolution Methods}
\usage{
  setGEDMethod(key, ...)

  removeGEDMethod(key, ...)
}
\arguments{
  \item{key}{accession key.}

  \item{...}{extra registry fields describing the method,
  or arguments passed to
  \code{\link[pkgmaker]{pkgreg_remove}}.}
}
\description{
  \code{setGEDMethod} register a deconvolution algorithm in
  the \pkg{CellMix} registry.

  \code{removeGEDMethod} removes a deconvlution algorithm
  from the registry. Note that it does not delete any NMF
  algorithm that might be related to it.
}

