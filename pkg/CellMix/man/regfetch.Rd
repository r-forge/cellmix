\name{regfetch}
\alias{regfetch}
\title{Advanced Fetch for Registry}
\usage{
  regfetch(regobj, key, msg = NULL, error = TRUE,
    exact = FALSE, multiple = FALSE, verbose = FALSE)
}
\arguments{
  \item{regobj}{registry object.}

  \item{key}{entry accession key.}

  \item{msg}{error prefix.}

  \item{error}{logical that indicates if an error should be
  thrown if the key is not found or if matches multiple
  entries (when \var{multiple=FALSE}).}

  \item{exact}{logical that indicates if the one looks for
  exact match(es).}

  \item{multiple}{logical that indicates if the key can
  return multiple matches, or \code{NA} if
  \code{error=FALSE}, or throw an error if
  \code{error=TRUE}.}

  \item{verbose}{logical that toggles verbosity.}
}
\description{
  Retrieve entries from a registry object.
}
\keyword{internal}

