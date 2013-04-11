\docType{methods}
\name{sapply,MarkerList-method}
\alias{sapply,MarkerList-method}
\title{Aplying Functions Along MarkerList Objects}
\usage{
  \S4method{sapply}{MarkerList}(X, FUN, ..., ASNEW = FALSE,
    ARG.NAMES = FALSE, simplify = TRUE, USE.NAMES = TRUE)
}
\arguments{
  \item{ASNEW}{logical that indicates if the returned
  object should have all data other than the list identical
  to \code{X}.}

  \item{ARG.NAMES}{logical that indicates if the names of
  each element should be passed to \code{FUN} as a second
  argument. In this case, the function is applie using
  \code{\link{mapply}}, with two vectorised arguments: the
  list itself and its names.}

  \item{simplify}{logical that indicates if the result
  should be simplified if possible. If \code{TRUE}, then
  the result is simplified as in
  \code{\link[base]{sapply}}. If \code{FALSE}, then the
  result is not simplified, and returned as a
  \code{\linkS4class{MarkerList}} object if possible, or as
  a list otherwise. If \code{NA}, then the result is not
  simplified, and always returned as a list.}

  \item{X}{ A list-like or vector-like R object.  }

  \item{FUN}{ See \code{?base::\link[base]{lapply}} for a
  description of these arguments.  }

  \item{...}{ See \code{?base::\link[base]{lapply}} for a
  description of these arguments.  }

  \item{USE.NAMES}{ See \code{?base::\link[base]{lapply}}
  for a description of these arguments.  }
}
\description{
  Applies a given function to a
  \code{\linkS4class{MarkerList}} object.
}
\keyword{methods}

