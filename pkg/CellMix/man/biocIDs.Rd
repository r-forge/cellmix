\name{is.probeid}
\alias{asGeneIdentifierType}
\alias{is.idtype}
\alias{is.probeid}
\alias{is.probetype}
\title{Utility function for Biological Identifiers}
\usage{
  is.probeid(x, ...)

  is.probetype(type)

  is.idtype(type)

  asGeneIdentifierType(x, limit = NULL, annotation = NULL,
    error = NullIdentifier(annotation))
}
\arguments{
  \item{x}{an R object}

  \item{...}{extra arguments passed to
  \code{\link{idtype}}.}

  \item{type}{an identifier type as returned by
  \code{\link{idtype}}.}

  \item{annotation}{annotation package to associate with
  the returned \code{\linkS4class{GeneIdentifierType}}
  object.}

  \item{error}{specifies what to do if the type of \code{x}
  cannot be determined

  For \code{asGeneIdentifierType}, all arguments in
  \code{...} are used in an internal call to
  \code{idtype}.}

  \item{limit}{specification for limiting which elements
  are used to detect the type of identifiers. If a single
  numeric, then only the first \code{limit} elements are
  used. Otherwise it must be a subsetting logical or
  numeric vector.}
}
\description{
  \code{is.probeid} tells if given IDs are probe IDs.

  \code{is.probetype} tells if given types are probe ID
  types.

  \code{is.idtype} tells if a given character vector
  contains valid types.

  \code{asGeneIdentifierType} is similar to
  \code{\link{idtype}}, but returns an object of class
  \code{\linkS4class{GeneIdentifierType}}, for use in
  \code{\linkS4class{GeneSet}} objects.
}

