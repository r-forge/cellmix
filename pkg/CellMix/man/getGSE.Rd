\name{getGSE}
\alias{getGSE}
\title{Downloading GSE Datasets from GEO}
\usage{
  getGSE(GEO = NULL, filename = NULL, destdir = GEDtmp(),
    ..., simplify = TRUE)
}
\arguments{
  \item{...}{extra parameters passed to
  \code{\link[GEOquery]{getGEO}}.}

  \item{simplify}{logical that indicates if the result for
  GSEs that contain only one dataset should be simplified
  and returned as an \code{ExpressionSet} object
  (\code{TRUE}) or as a list.}

  \item{GEO}{A character string representing a GEO object
  for download and parsing.  (eg.,
  'GDS505','GSE2','GSM2','GPL96')}

  \item{filename}{The filename of a previously downloaded
  GEO SOFT format file or its gzipped representation (in
  which case the filename must end in .gz).  Either one of
  GEO or filename may be specified, not both.  GEO series
  matrix files are also handled.  Note that since a single
  file is being parsed, the return value is not a list of
  esets, but a single eset when GSE matrix files are
  parsed.}

  \item{destdir}{The destination directory for any
  downloads.  Defaults to the architecture-dependent
  tempdir.  You may want to specify a different directory
  if you want to save the file for later use.  Doing so is
  a good idea if you have a slow connection, as some of the
  GEO files are HUGE!}
}
\description{
  This function acts as a wrapper for
  \code{\link[GEOquery]{getGEO}}, to download the matrix
  files of gene expression datasets from the \code{GEO}
  database. It fixes two issues:
}
\details{
  \itemize{ \item if the user is behind a proxy that
  modifies the content of GEO's dataset index plain text
  pages into HTML pages, then plain \code{getGEO} breaks
  with an error like:

  \code{Error in scan(file, what, nmax, sep, dec, quote,
  skip, nlines, na.strings, : line 1 did not have 6
  elements} \item plain \code{getGEO} do not honour the
  destination directory for GPL files. }

  Moreover, this function installs GEOquery if necessary.

  This function is a wrapper to
  \code{\link[GEOquery]{getGEO}}, which unlist the result
  list in the common case of a single dataset, and installs
  the \pkg{GEOquery} package if necessary.
}
\examples{
\dontrun{
dir.create('testGSE')
getGSE('GSE19830', destdir='testGSE')
}
}

