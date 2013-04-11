\docType{methods}
\name{convertIDs}
\alias{convertIDs}
\alias{convertIDs,ANY,ANY,ANY-method}
\alias{convertIDs,ANY,ANY,NullIdentifier-method}
\alias{convertIDs,ANY,list,missing-method}
\alias{convertIDs,character,GeneIdentifierType,GeneIdentifierType-method}
\alias{convertIDs,ExpressionSet,GeneIdentifierType,GeneIdentifierType-method}
\alias{convertIDs,list,GeneIdentifierType,GeneIdentifierType-method}
\alias{convertIDs,MarkerList,GeneIdentifierType,GeneIdentifierType-method}
\alias{convertIDs,matrix,GeneIdentifierType,GeneIdentifierType-method}
\alias{convertIDs-methods}
\alias{mapIdentifiers,list,GeneIdentifierType,GeneIdentifierType-method}
\title{Converting Gene or Probeset IDs}
\usage{
  convertIDs(object, to, from, ...)

  \S4method{mapIdentifiers}{list,GeneIdentifierType,GeneIdentifierType}(what,
    to, from, ..., verbose = FALSE)

  \S4method{convertIDs}{MarkerList,GeneIdentifierType,GeneIdentifierType}(object,
    to, from, verbose = FALSE, nodups = NULL,
    unlist = TRUE, ...)

  \S4method{convertIDs}{character,GeneIdentifierType,GeneIdentifierType}(object,
    to, from, method = "auto", unlist = TRUE, ...)

  \S4method{convertIDs}{matrix,GeneIdentifierType,GeneIdentifierType}(object,
    to, from, ..., unlist = TRUE, rm.duplicates = NULL)

  \S4method{convertIDs}{ExpressionSet,GeneIdentifierType,GeneIdentifierType}(object,
    to, from, ..., unlist = TRUE, rm.duplicates = NULL)
}
\arguments{
  \item{object}{object whose identifiers are converted.}

  \item{to}{specification of the type of identifiers to
  convert to.}

  \item{from}{specification of the type of identifiers of
  \code{object}. This is only neeeded when the source type
  cannot be inferred from \code{object} itself.}

  \item{...}{extra arguments to allow extension, which are
  passed down to the workhorse method
  \code{convertIDs,character,GeneIdentifierType,GeneIdentifierType}.
  See each method's description for more details.}

  \item{what}{identifiers to map}

  \item{verbose}{a logical or integer that sets the
  vverbosity level.}

  \item{nodups}{specifies if marker identifiers that are
  duplicated across cell types in the result should be
  removed (\code{TRUE}) or not (\code{FALSE}). If
  \code{NULL}, then duplicates are removed only if there
  were no duplicates in the source object.}

  \item{method}{mapping method, passed to
  \code{\link{mapIDs}}, that indicates how to carry the
  mapping between the original and final identifier type.}

  \item{unlist}{logical that indicates if the result should
  be flatten, i.e. turned into a vector rather than a list
  -- using \code{\link{unlist2}}. In this case, the
  vector's name then correspond to the source identifiers.}

  \item{rm.duplicates}{logical that how duplicated should
  be treated. \code{rm.duplicates=FALSE} does not allow any
  duplicated match and throws an error if any is present.
  If \code{TRUE} or \code{NULL} duplicates only the first
  match is kept, but a warning is thrown only when
  \code{NULL}.}
}
\description{
  The \pkg{CellMix} package implements methods for
  converting markers identifiers cross-platform and
  cross-organism, to facilitate deconvolution analysis to
  be carried out using multiple independent sources of
  data, e.g., use a marker gene list obtained from data on
  one platform to deconvolve gene expression data generated
  on another platform.
}
\details{
  The function \code{convertIDs} provides the main
  interface to convert genes/probeset ids into IDs
  compatible with another given data. It is typically
  useful to convert built-in marker gene lists (see
  \code{link{cellMarkers}}).

  The identifier conversion functions and methods defined
  in the \pkg{CellMix} package can be seen as extending the
  existing framework defined in the \pkg{GSEABase} package,
  with the generic \code{\link[GSEABase]{mapIdentifiers}}.
}
\section{Methods}{
  \describe{

  \item{convertIDs}{\code{signature(object = "list", to =
  "GeneIdentifierType", from = "GeneIdentifierType")}:
  Apply the conversion to each element of the list. }

  \item{convertIDs}{\code{signature(object = "MarkerList",
  to = "GeneIdentifierType", from = "GeneIdentifierType")}:
  Convert IDs from a MarkerList object.

  In this case, argument \code{unlist} indicates if the
  result should be a simple list containing the mapping (a
  list) for each cell type or a
  \code{\linkS4class{MarkerList}} object (default). }

  \item{convertIDs}{\code{signature(object = "character",
  to = "GeneIdentifierType", from = "GeneIdentifierType")}:
  This is the workhorse method that is eventually called by
  all other \code{convertIDs} methods. The actual
  conversions are perforemd by \code{\link{mapIDs}}, to
  which are passed all arguments in \code{...}, in
  particular, arguments \code{verbose} and \code{method}. }

  \item{convertIDs}{\code{signature(object = "matrix", to =
  "GeneIdentifierType", from = "GeneIdentifierType")}:
  Convert the row names of a matrix into other identifiers.

  In this case, argument \code{unlist} indicates if the
  converted ids should be used to subset the original
  \code{matrix} object, or returned directly returned as a
  list. }

  \item{convertIDs}{\code{signature(object =
  "ExpressionSet", to = "GeneIdentifierType", from =
  "GeneIdentifierType")}: Convert the feature names of an
  ExpressionSet into other identifiers.

  In this case, argument \code{unlist} indicates if the
  converted ids should be used to subset the original
  \code{ExpressionSet} object, or returned directly
  returned as a list. }

  \item{convertIDs}{\code{signature(object = "ANY", to =
  "ANY", from = "NullIdentifier")}: Convert identifiers,
  inferring the type of origin from the object itself, but
  keep the annotation specification embedded in
  \code{from}. }

  \item{convertIDs}{\code{signature(object = "ANY", to =
  "ANY", from = "ANY")}: Convert identifiers, inferring the
  type from the specifications in \code{to} and
  \code{from}, eg., \code{to='ENTREZID'}, or
  \code{'UNIGENE'}. If not specified in either \code{to} or
  \code{from}, the annotation is taken from \code{object}.
  If \code{from} is missing, the source type is infered
  from \code{object} itself. }

  \item{convertIDs}{\code{signature(object = "ANY", to =
  "list", from = "missing")}: Convert identifiers using a
  given map or list of maps. }

  \item{mapIdentifiers}{\code{signature(what = "list", to =
  "GeneIdentifierType", from = "GeneIdentifierType")}:
  Applies \code{mapIdentifier} to each element in a list.

  All arguments in \code{...} are passed to the subsequent
  calls to \code{\link[GSEABase]{mapIdentifiers}}. }

  }
}
\examples{
# load a marker list from the registry
m <- MarkerList('IRIS')
summary(m)
head(m[[1]])

# convert Entrez gene ids to Affy probeset ids chip hgu133b
m2 <- convertIDs(m, 'hgu133b.db', verbose=2)
summary(m2)
#----------------------------------------------
# 1. Conversion from biological IDs
#----------------------------------------------
# For this kind of IDs, a source annotation package can often be inferred
# from the ID type, using regular expression patterns (e.g. "^ENS[0-9]+$"
# identifies Ensembl gene IDs)

ids <- c("Hs.1", "Hs.2", "Hs.3")
# get Entrez gene IDs (based on annotation from the org.Hs.gene.eg package)
convertIDs(ids, 'ENTREZID', 'org.Hs.eg.db', verbose=TRUE)

# map to other IDs
convertIDs(ids, 'REFSEQ')
convertIDs(ids, 'ENSEMBL')
# convert across ogranism
convertIDs(ids, 'rat2302.db')
# get Affy probeset IDs for chip hgu133a
affy <- convertIDs(ids, 'hgu133a.db')

# assume we have a vector of IDs, e.g. Entrez gene ids
id <- c("673", "725", "10115")
# get associated probesets on chip hgu133a
convertIDs(id, 'hgu133a.db')
# get all associated probesets on chip hgu133a
convertIDs(id, 'hgu133a.db', method='all')
# same as a vector with duplicated names
convertIDs(id, 'hgu133a.db', method='all', unlist=FALSE)
# specification using ProbeAnnDbBimap objects
library(hgu133b.db)
convertIDs(id, 'hgu133a.db', hgu133bENTREZID, verbose=2)

#----------------------------------------------
# 2. Conversion from probeset IDs
#----------------------------------------------
# For this kind of IDs, a source annotation package is required, because it
# cannot be easily inferred from the ID type.

# get Affy probeset IDs for chip hgu133b from ids for hgu133b
convertIDs(affy, 'hgu133a.db', 'hgu133b.db')
# across organism
convertIDs(affy, 'hgu133a.db', 'rat2302.db')
}
\keyword{methods}

