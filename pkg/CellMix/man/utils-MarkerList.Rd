\docType{methods}
\name{dropvalues}
\alias{drop,MarkerList-method}
\alias{dropvalues}
\alias{hasDuplicated}
\alias{mltype}
\alias{nmark}
\alias{nmark,ANY-method}
\alias{nmark,list-method}
\alias{nmark,MarkerList-method}
\alias{nmark-methods}
\alias{reverse}
\alias{reverse,MarkerList-method}
\alias{reverse-methods}
\alias{rmDuplicated}
\title{Utility Functions for MarkerList Objects}
\usage{
  dropvalues(object, ..., int.ok = TRUE)

  mltype(object,
    type = c("character", "numeric", "integer", "logical"))

  rmDuplicated(object, ...)

  hasDuplicated(object, which = FALSE)

  reverse(object, ...)

  \S4method{reverse}{MarkerList}(object, all = FALSE)

  nmark(object, ...)

  \S4method{nmark}{list}(object, each = FALSE)

  \S4method{drop}{MarkerList}(x)
}
\arguments{
  \item{object}{an R object, typically a list that contains
  marker data structured in a similar way as in
  \code{MarkerList} objects.}

  \item{...}{extra argument to allow extension}

  \item{int.ok}{logical that indicates if integer marker
  names are valid, or if one should use their names
  instead.}

  \item{type}{type to compare with.}

  \item{which}{logical that indicates if a list of
  duplicated identifiers should be returned. If
  \code{which=TRUE} and no duplicates exist then
  \code{NULL} is returned.}

  \item{all}{a logical that indicates if the list itself,
  i.e., the sets should be reversed as well.}

  \item{each}{a logical indicating if the number of marker
  in each set is to be returned or only the total number
  (default).}

  \item{x}{an array (including a matrix).}
}
\description{
  \code{dropvalues} Drops the values associated with each
  marker (e.g. tissue specificity score), and returns an
  object of the same type as the input object. It is a
  shorcut for \code{\link{marknames}(object,
  unlist=FALSE)}.

  \code{mltype} returns the type of marker list or compare
  it with a given type: \itemize{ \item character: elements
  are marker names/ids; \item numeric: elements are numeric
  vectors (e.g., scores), named with the marker names;
  \item index: elements are integer vectors corresponding
  to indexes in some reference data; \item logical:
  elements are logical vectors, named with marker names,
  which result from logical operations on a
  \code{MarkerList} object. }

  \code{rmDuplicated} remove marker identifiers that are
  duplicated either within or between sets. Arguments in
  \code{...} are not used.

  \code{hasDuplicated} checks for duplicated identifiers
  across sets.

  \code{reverse} reverse is an S4 method to reverse an
  object. For \code{MarkerList} object, it reverse the
  order of markers within each set seprately, and
  optionally the order of the set in the list.

  \code{nmark} returns the number of markers for each
  cell-type. It accepts any object as an input, but is
  generally called on an \code{MarkerList} object or an
  object with an attached \code{MarkerList} object (see
  methods \code{getMarkers} and \code{attachMarkers}).
}
\section{Methods}{
  \describe{

  \item{drop}{\code{signature(x = "MarkerList")}:
  \code{drop} drops empty sets of markers from a
  \code{MarkerList} object, as well as markers with
  \code{FALSE} values if the object is a logical marker
  list. It returns the reduced \code{MarkerList} object. }

  \item{nmark}{\code{signature(object = "list")}: The
  method for lists, including \code{MarkerList} objects
  returns the total number of markers, or optionaly its
  breakdown per set of markers. }

  \item{nmark}{\code{signature(object = "ANY")}: The
  default method tries to extract embedded marker data from
  \code{object} using \code{\link{getMarkers}}, and returns
  the result of \code{nmark} applied to these data if
  present or 0L otherwise. }

  }
}
\examples{
m <- gMarkerList(3, 4)
geneIds(m)
geneIds( reverse(m) )
geneIds( reverse(m, all=TRUE) )
m <- MarkerList( list(1:4, 10:15, integer(), 100) )
geneIds(m)
geneIds(drop(m))
}
\keyword{methods}

