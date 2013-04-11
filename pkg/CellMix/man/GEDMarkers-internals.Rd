\name{setMarkerList}
\alias{removeMarkerList}
\alias{setMarkerList}
\title{Register CellMix MarkerLists}
\usage{
  setMarkerList(key, ...)

  removeMarkerList(key, ...)
}
\arguments{
  \item{key}{accession key.}

  \item{...}{extra registry fields describing the marker
  list, or arguments passed to
  \code{\link[pkgmaker]{pkgreg_remove}}.}
}
\description{
  \code{setMarkerList} register a marker list in the
  \pkg{CellMix} registry.

  \code{removeMarkerList} removes a marker list from the
  registry.
}

