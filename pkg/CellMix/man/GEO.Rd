\name{getAndParseGSEMatrices}
\alias{getAndParseGSEMatrices}
\alias{parseGSEMatrix}
\title{Custom GEO Functions}
\usage{
  getAndParseGSEMatrices(GEO, destdir, AnnotGPL)

  parseGSEMatrix(fname, AnnotGPL = FALSE,
    destdir = GEOtempdir())
}
\arguments{
  \item{GEO}{geo accession key}

  \item{destdir}{destination directory}

  \item{AnnotGPL}{logical that indicates if GPL data should
  be downloaded to annotate the features.}

  \item{fname}{filename of the GSE matrix file}
}
\description{
  \code{getAndParseGSEMatrices} is an internal modified
  version of \code{GEOquery::getAndParseGSEMatrices}, which
  works even when the use is behind a proxy which turns
  GEO's dataset index plain text pages into HTML pages.
  Moreover it correctly passes the destination directory
  argument (\code{destdir}) to the subsequent calls, so
  that GPL description files are also downloaded there.

  \code{parseGSEMatrix} fixes the issue that the
  destination directory for GPL files is not passed down to
  other getGEO calls.
}
\keyword{internal}

