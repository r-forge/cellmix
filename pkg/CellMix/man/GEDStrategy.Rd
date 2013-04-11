\docType{methods}
\name{GEDStrategy}
\alias{GEDStrategy}
\alias{GEDStrategy,character-method}
\alias{GEDStrategy,function-method}
\alias{GEDStrategy,GEDalgorithm_entry-method}
\alias{GEDStrategy,GEDStrategy-method}
\alias{GEDStrategy-methods}
\alias{GEDStrategy,missing-method}
\title{Factory Methods for GEDStrategy Objects}
\usage{
  GEDStrategy(object, ...)

  \S4method{GEDStrategy}{GEDStrategy}(object, key, ...)

  \S4method{GEDStrategy}{function}(object, ...,
    name = NULL)
}
\arguments{
  \item{object}{object from which of a \code{GEDalgorithm}
  is created}

  \item{...}{extra arguments to allow extension}

  \item{key}{if not missing, it must be a character string,
  which is used as unique key for the new GED algorithm
  object. In particular, it cannot be a key already used by
  one of the registered algorithms, except if \code{object}
  and \code{key} are the only arguments and if the key
  matches the key of \code{object} itself, which simply
  returns the object unchanged -- and is somehow not a very
  useful call. If missing, a unique key is generated using
  a digest hash of the result from
  \code{\link{tempfile}()}.}

  \item{name}{descriptive name of the algorithm. A
  temporary name is generated if missing.}
}
\description{
  Factory Methods for GEDStrategy Objects
}
\section{Methods}{
  \describe{

  \item{GEDStrategy}{\code{signature(object = "missing")}:
  This method serves as a standard constructor and is
  equivalent to \code{new('GEDalgorithm', ...)}, meaning
  that all arguments must be named and correspond to slots
  in \code{\linkS4class{GEDStrategy}}. }

  \item{GEDStrategy}{\code{signature(object =
  "GEDalgorithm_entry")}: Method to create a
  \code{GEDalgorithm} object from its entry in the internal
  registry. }

  \item{GEDStrategy}{\code{signature(object =
  "GEDStrategy")}: Constructor-Copy method, which creates
  an \code{GEDalgorithm} object from a template and ensure
  that the object has a key that does not conflict with
  registered algorithms. }

  \item{GEDStrategy}{\code{signature(object =
  "character")}: This method creates a \code{GEDalgorithm}
  object from an access key, which allows to make objects
  using any registered algorithm as a template.

  The access key must be valid, i.e. that it must
  correspond exactly to a registered algorithm. If other
  arguments are passed, these must correspond to slots in
  \code{\linkS4class{GEDStrategy}}, and will overwrite the
  values defined by the registered algorithm. }

  \item{GEDStrategy}{\code{signature(object = "function")}:
  This method creates a \code{GEDalgorithm} object from a
  plain R function. }

  }
}
\keyword{methods}

