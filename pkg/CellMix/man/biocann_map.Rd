\name{biocann_map}
\alias{biocann_map}
\title{Retrieving ID Conversion Maps}
\usage{
  biocann_map(from, to, keys = NULL, cache = TRUE,
    link = "ENTREZID", lib = NULL, verbose = FALSE)
}
\arguments{
  \item{from}{the source tyep of gene identifier specified
  as a \code{\linkS4class{GeneIdentifierType}} object, or a
  character string, e.g. \code{'ENTREZID'}, if \code{to}
  already contains annotation data.}

  \item{to}{the source tyep of gene identifier specified as
  a \code{\linkS4class{GeneIdentifierType}} object, or a
  character string, e.g. \code{'UNIGENE'}, if \code{from}
  already contains annotation data.}

  \item{keys}{character vector of gene identifiers to map.
  This is only used when performing cross-species mapping.}

  \item{cache}{logical that indicates if -- lengthy --
  cross-species mappings should be cached, so that
  subsequent calls with the same set of keys returns much
  faster.}

  \item{link}{the gene identifier type to use to link
  between platforms. We \strong{strongly} recommend to
  leave this argument as default \code{'ENTREZID'}.}

  \item{lib}{path to the library (ie. directory) where to
  first look and install any required annotation packages.
  This library is temporarly prepended to the .libPath, so
  that user and R standard library locations are still
  looked up.}

  \item{verbose}{logical or integer that indicates the
  level of verbosity.}
}
\description{
  Computes a complete list of maps to apply in order to
  convert from one type of gene identifier to another.
}
\details{
  This function is an extension/enhancement of the function
  \code{GSEABase:::.mapIdentifiers_selectMaps}, which cope
  with the cases of cross-platform and cross-species
  conversions.

  In particular, when called in an interactive session,
  this function will notifies the user about missing
  required annotation packages, and offers to install them
  at runtime.
}
\examples{
\dontrun{
# intra-platform
m <- biocann_map('ENTREZID', AnnotationIdentifier('hgu133a.db'))
# cross platform
m <- biocann_map(UnigeneIdentifier('hgu133a.db'), EntrezIdentifier('hgu133b.db'))
# cross platform, cross-species
m <- biocann_map(EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('mogene10stprobeset.db'), keys=as.character(1:1000))
}
}

