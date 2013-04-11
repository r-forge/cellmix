\docType{methods}
\name{geneValues}
\alias{annotation}
\alias{annotation<-}
\alias{annotation,GEDdata_entry-method}
\alias{annotation<-,MarkerList,character-method}
\alias{annotation,MarkerList-method}
\alias{annotation<-,MarkerList,NULL-method}
\alias{annotation<--methods}
\alias{annotation-methods}
\alias{details}
\alias{details,MarkerList-method}
\alias{details-methods}
\alias{geneIds}
\alias{geneIds<-}
\alias{geneIds<-,MarkerList,list-method}
\alias{geneIds,MarkerList-method}
\alias{geneIds<--methods}
\alias{geneIds-methods}
\alias{geneIdType}
\alias{geneIdType<-}
\alias{geneIdType<-,MarkerList,character-method}
\alias{geneIdType<-,MarkerList,GeneIdentifierType-method}
\alias{geneIdType,MarkerList-method}
\alias{geneIdType<-,MarkerList,NULL-method}
\alias{geneIdType<--methods}
\alias{geneIdType-methods}
\alias{geneValues}
\alias{geneValues<-}
\alias{geneValues,MarkerList-method}
\alias{geneValues<--methods}
\alias{geneValues-methods}
\alias{hasValues}
\alias{marknames}
\alias{marknames,list-method}
\alias{marknames-methods}
\alias{marknames,vector-method}
\alias{nmf}
\alias{nmf,MatrixData,MarkerList,ANY-method}
\alias{nmf-methods}
\alias{show,MarkerList-method}
\title{Accessing Data in Marker Lists}
\usage{
  geneValues(object, ...)

  geneValues(object, ...) <- value

  geneIds(object, ...)

  geneIds(object) <- value

  geneIdType(object)

  geneIdType(object, verbose = FALSE) <- value

  annotation(object, ...)

  annotation(object, ...) <- value

  details(object)

  nmf(x, rank, method, ...)

  \S4method{show}{MarkerList}(object)

  \S4method{details}{MarkerList}(object)

  \S4method{geneValues}{MarkerList}(object, unlist = FALSE)

  hasValues(object, ...)

  marknames(object, ...)

  \S4method{marknames}{vector}(object, simplify = TRUE)

  \S4method{marknames}{list}(object, unlist = TRUE, ...)
}
\arguments{
  \item{object}{an object from/which to get/set data}

  \item{value}{replacement value}

  \item{...}{extra arguments to allow extension. See each
  method's description for more details.}

  \item{verbose}{logical not used.}

  \item{unlist}{a logical that indicates if the values
  should be returned as a list (\code{FALSE}) or as a
  vector.}

  \item{simplify}{logical that is use when for logical
  marker lists only, and indicates if all marker names
  should be returned (\code{simplify=FALSE}), or only the
  ones associated with a \code{TRUE} value.}

  \item{x}{target data to fit, i.e. a matrix-like object}

  \item{rank}{specification of the factorization rank. It
  is usually a single numeric value, but other type of
  values are possible (e.g. matrix), for which specific
  methods are implemented. See for example methods
  \code{nmf,matrix,matrix,ANY}.

  If \code{rank} is a numeric vector with more than one
  element, e.g. a range of ranks, then \code{\link{nmf}}
  performs the estimation procedure described in
  \code{\link{nmfEstimateRank}}.}

  \item{method}{specification of the NMF algorithm. The
  most common way of specifying the algorithm is to pass
  the access key (i.e. a character string) of an algorithm
  stored in the package's dedicated registry, but methods
  exists that handle other types of values, such as
  \code{function} or \code{list} object. See their
  descriptions in section \emph{Methods}.

  If \code{method} is missing the algorithm to use is
  obtained from the option
  \code{nmf.getOption('default.algorithm')}, unless it can
  be infer from the type of NMF model to fit, if this later
  is available from other arguments. Factory fresh default
  value is \sQuote{brunet}, which corresponds to the
  standard NMF algorithm from \cite{Brunet2004} (see
  section \emph{Algorithms}).

  Cases where the algorithm is inferred from the call are
  when an NMF model is passed in arguments \code{rank} or
  \code{seed} (see description for
  \code{nmf,matrix,numeric,NULL} in section
  \emph{Methods}).}
}
\value{
  a character vector, usually of length 1 but possibly
  longer.
}
\description{
  The S4 generics \code{geneValues} and \code{geneValues<-}
  respectively return and set the values associated with
  each gene in a gene set.

  \code{geneIds} returns all gene identifiers in a
  \code{MarkerList} object as a standard list.

  \code{annotation} and \code{annotation<-} respectively
  extract and set the annotation string from a
  \code{\link{MarkerList}} object, i.e. the name of
  Bioconductor annotation package(s) relevant to convert
  marker identifiers into other IDs.

  \code{hasValues} tells if a marker list contains numeric
  values, e.g. specificity scores.

  \code{marknames} is similar to \code{geneIds} but by
  default returns the all marker names unlisted, i.e. as a
  character vector.
}
\section{Methods}{
  \describe{

  \item{annotation}{\code{signature(object =
  "MarkerList")}: Extract name of the annotation package
  from a marker list. }

  \item{annotation}{\code{signature(object =
  "GEDdata_entry")}: \code{annotation} gets the name of the
  annotation package(s) relevant for the data. }

  \item{annotation<-}{\code{signature(object =
  "MarkerList", value = "character")}: Sets the name of the
  annotation package associated to a marker list. }

  \item{annotation<-}{\code{signature(object =
  "MarkerList", value = "NULL")}: Reset the annotation
  package associated to a marker list }

  \item{details}{\code{signature(object = "MarkerList")}:
  Shows more details than the regular \code{show} method. }

  \item{geneIdType}{\code{signature(object =
  "MarkerList")}: Returns the type of gene identifier used
  in a marker list. }

  \item{geneIdType<-}{\code{signature(object =
  "MarkerList", value = "NULL")}: Sets the type of gene
  identifier used in a marker list to
  \code{\link{NullIdentifier}()}. }

  \item{geneIdType<-}{\code{signature(object =
  "MarkerList", value = "character")}: Sets the type of
  gene identifier used in a marker list to a given value. }

  \item{geneIdType<-}{\code{signature(object =
  "MarkerList", value = "GeneIdentifierType")}: Sets the
  type of gene identifier associated with a marker list to
  a given \code{\linkS4class{GeneIdentifierType}} object. }

  \item{geneValues}{\code{signature(object =
  "MarkerList")}: Returns the values embedded in a
  \code{MarkerList} object. }

  \item{nmf}{\code{signature(x = "MatrixData", rank =
  "MarkerList", method = "ANY")}: Calls an NMF algorithm
  using a MarkerList object, whose length defines the
  factorization rank. The actual MarkerList object is
  passed down in argument \code{data}. }

  }
}
\examples{
# get marker Ids
m <- rMarkerList(3, names=TRUE)
geneIds(m)
# get attached values
m <- rMarkerList(3, values=TRUE)
geneValues(m)
marknames( rMarkerList(3) )
marknames( rMarkerList(3, names=TRUE) )
marknames( rMarkerList(3, values=TRUE) )
}
\keyword{methods}

