\docType{methods}
\name{idtype}
\alias{idtype}
\alias{idtype,ChipDb-method}
\alias{idtype,ExpressionSet-method}
\alias{idtype,GeneIdentifierType-method}
\alias{idtype,list-method}
\alias{idtype,MarkerList-method}
\alias{idtype,matrix-method}
\alias{idtype-methods}
\alias{idtype,missing-method}
\alias{idtype,NMF-method}
\alias{idtype,NULL-method}
\alias{idtype,ProbeAnnDbBimap-method}
\alias{idtype,vector-method}
\title{Identifying Gene or Probe ID Type}
\usage{
  idtype(object, ...)

  \S4method{idtype}{missing}(object, def = FALSE)

  \S4method{idtype}{ProbeAnnDbBimap}(object, limit = 500L,
    ...)

  \S4method{idtype}{ChipDb}(object, limit = 500L, ...)

  \S4method{idtype}{MarkerList}(object, each = FALSE, ...)

  \S4method{idtype}{vector}(object, each = FALSE,
    limit = NULL, no.match = "")
}
\arguments{
  \item{object}{an R object that contains the gene
  identifiers whose type is to be determined.}

  \item{...}{extra argument to allow extension, generally
  passed down to \code{idtype,character-method}. See each
  method's description for more details.}

  \item{def}{a logical or a subsetting vector, used when
  \code{object} is missing, which indicates that the result
  should contain the definition of the matching
  pattern/function of each type, or which type's deifnition
  should be included in the result list.}

  \item{each}{logical indicating whether the type of each
  element should be returned (\code{TRUE}) or only the type
  of the vector as a whole (default).}

  \item{limit}{specification for limiting which elements
  are used to detect the type of identifiers. If a single
  numeric, then only the first \code{limit} elements are
  used. Otherwise it must be a subsetting logical or
  numeric vector.}

  \item{no.match}{character string that specifies the
  string to use when the type cannot be determined.

  The IDs can be either: \itemize{ \item{probe IDs (e.g.
  123456_at or ILMN_123456 for Affymetrix or Illumina chips
  respectively), the type starts with a dot \code{'.'},
  allowing the subsequent handling of such IDs as a group.}
  \item{other biological ID types, the result are character
  strings such as those used as attributes in Bioconductor
  annotation packages (e.g. \code{"ENTREZID"} or
  \code{"ENSEMBL"})} \item{Names of annotation packages
  e.g. \code{"hgu133plus2.db"}.} }

  This function is able to identify the following ID types
  using regular expression patterns or dedicated function:
  \itemize{ \item ENSEMBL = "^ENSG[0-9]+$" \item
  ENSEMBLTRANS = "^ENST[0-9]+$" \item ENSEMBLPROT =
  "^ENSP[0-9]+$" \item ENTREZID = "^[0-9]+$" \item IMAGE =
  "^IMAGE:[0-9]+$" \item GOID = "^GO:[0-9]+$" \item PFAM =
  "^PF[0-9]+$" \item REFSEQ = "^N[MP]_[0-9]+$" \item ENZYME
  = "^[0-9]+(\\.(([0-9]+)|-)+){3}$" \item MAP =
  "^[0-9XY]+((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9XY]+)?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$"
  \item GENEBANK (Nucleotide) = "^[A-Z][0-9]{5}$" |
  "^[A-Z]{2}[0-9]{6}$" \item GENEBANK (Protein) =
  "^[A-Z]{3}[0-9]{5}$" \item GENEBANK (WGS) =
  "^[A-Z]{4}[0-9]{8}[0-9]?[0-9]?$" \item GENEBANK (MGA) =
  "^[A-Z]{5}[0-9]{7}$" \item GENENAME = " " \item
  .Affymetrix =
  "(^AFFX-)|(^[0-9]+_([abfgilrsx]_)?([as]t)|(i))$" \item
  .Illumina = "^ILMN_[0-9]+$" \item .Agilent =
  "^A_[0-9]+_P[0-9]+$" \item .nuID = use the function
  \code{\link{nuIDdecode}} to try converting the ids into
  nucleotide sequences. Identification is positive if no
  error is thrown during the conversion. }}
}
\value{
  a single character string (possibly empty) if
  \code{each=FALSE} (default) or a character vector of the
  same "length" as \code{object} otherwise.
}
\description{
  The S4 generic \code{idtype} automatically determine the
  type of gene/feature identifiers stored in objects, based
  on a combination of regular expression patterns and test
  functions.
}
\details{
  It uses a heuristic based on a set of regular expressions
  and functions that uniquely match most common types of
  identifiers, such as Unigene, entrez gene, Affymetrix
  probe ids, Illumina probe ids, etc..
}
\section{Methods}{
  \describe{

  \item{idtype}{\code{signature(object = "missing")}:
  Method for when \code{idtype} is called with its first
  argument missing, in which case it returns all or a
  subset of the known type names as a character vector, or
  optionally as a list that contains their definition, i.e.
  a regular expression or a matching function. }

  \item{idtype}{\code{signature(object = "matrix")}:
  Detects the type of identifiers used in the row names of
  a matrix. }

  \item{idtype}{\code{signature(object = "ExpressionSet")}:
  Detects the type of identifiers used in the feature names
  of an \code{ExpressionSet} object. }

  \item{idtype}{\code{signature(object = "NMF")}: Detects
  the type of identifiers used in the rownames of the basis
  matrix of an \code{NMF} model. }

  \item{idtype}{\code{signature(object =
  "ProbeAnnDbBimap")}: Detects the type of the primary
  identifiers of a probe annotation bimap object.

  To speedup the identification, only the first 500 probes
  are used by default, since the IDs are very likely to
  have been curated and to be of the same type. This can be
  changed using argument \code{limit}. }

  \item{idtype}{\code{signature(object = "ChipDb")}:
  Detects the type of the identifiers of a chip annotation
  object.

  To speedup the identification, only the first 500 probes
  are used by default, since the IDs are very likely to
  have been curated and to be of the same type. This can be
  changed using argument \code{limit}. }

  \item{idtype}{\code{signature(object =
  "GeneIdentifierType")}: Returns the type of identifier
  defined by a \code{\link{GeneIdentifierType}} object.
  Note that this methods is a bit special in the sense that
  it will return the string \dQuote{ANNOTATION} for
  annotation based identifiers, but will not tell which
  platform it is relative to. This is different to what
  \code{idtype} would do if applied to the primary
  identifiers of the corresponding annotation package. }

  \item{idtype}{\code{signature(object = "list")}: Detects
  the type of all elements in a list, but provides the
  option of detecting the type of each element separately.
  }

  \item{idtype}{\code{signature(object = "NULL")}: Dummy
  method -- defined for convenience -- that returns
  \code{''} }

  \item{idtype}{\code{signature(object = "vector")}: This
  is the workhorse method that determine the type of ids
  contained in a character vector. }

  }
}
\examples{
# all known types
idtype()
# with their definitions
idtype(def=TRUE)
idtype(def='ENTREZID')
idtype(def=c('ENTREZID', 'ENSEMBLTRANS'))
# from GeneIdentifierType objects
idtype(NullIdentifier())
idtype(AnnotationIdentifier('hgu133a.db'))
# but
\dontrun{
	library(hgu133a.db)
	idtype(hgu133a.db)
}
idtype("12345_at")
idtype(c("12345_at", "23232_at", "555_x_at"))
# mixed types
ids <- c("12345_at", "23232_at", "Hs.1213")
idtype(ids) # not detected
idtype(ids, each=TRUE)
}
\keyword{methods}

