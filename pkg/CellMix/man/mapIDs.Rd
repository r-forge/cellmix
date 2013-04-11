\name{mapIDs}
\alias{mapIDs}
\title{Mapping Gene Identifiers}
\usage{
  mapIDs(keys, from, to,
    method = c("auto", "mauto", "all", "firstN", "1:1", "1:*", "affy"),
    mapper = biocann_map, verbose = FALSE, ...)
}
\arguments{
  \item{keys}{character vector of identifiers to map.}

  \item{method}{a specification of the filtering strategy.
  See \code{\link{idFilter}}.}

  \item{mapper}{mapper function, responsible to generate
  the actual mappings between identifiers. See details.}

  \item{from}{specification of the type of identifiers of
  \code{object}. This is only neeeded when the source type
  cannot be inferred from \code{object} itself.}

  \item{to}{specification of the type of identifiers to
  convert to.}

  \item{verbose}{a logical or integer that sets the
  vverbosity level.}

  \item{...}{extra arguments to allow extension, which are
  passed down to the workhorse method
  \code{convertIDs,character,GeneIdentifierType,GeneIdentifierType}.
  See each method's description for more details.}
}
\description{
  This function implements a generic mapping workflow that
  enables mapping gene identifers between different types,
  given a mapper function. The default mapper,
  \code{\link{biocann_map}}, enables mapping identifiers
  within or across platform, as well as across species.
}
\details{
  The mapper function passed in argument \code{mapper} is
  responsible for providing a sequence of map(s) that are
  sequentially applied, starting with the source gene
  identifiers in \code{keys}. It must at least have the
  following 4 arguments:

  \describe{ \item{from}{source gene identifier type, which
  will be passed as a \code{\link{GeneIdentifierType}})
  object.} \item{to}{source gene identifier type, which
  will be passed as a \code{\link{GeneIdentifierType}})
  object.} \item{keys}{the query keys to map, which may be
  used to help generate maps specific to the query, e.g.,
  for generating cross-species mappings. They may also be
  ignored, if already compiled maps exist, e.g., maps in
  standard annotation packages, as the returned map do not
  need to be limited to the query keys.} \item{verbose}{a
  logical or integer that toggles verbose messages at
  different levels.} }
}
\examples{
# some ENTREZ ids
ez <- as.character(1:1000)

# convert into probe ids for Affymetrix Chip hgu133a
m <- mapIDs(ez, EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('hgu133b.db'), verbose=2)

# keep primary affy probes only
m <- mapIDs(ez, EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('hgu133b.db'), method='affy', verbose=2)

# same but only keep 1:1 mapping, using a composed filtering strategy
m <- mapIDs(ez, EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('hgu133b.db'), method=c('affy', '1:1'), verbose=2)
}

