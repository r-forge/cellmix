\name{gedRequired}
\alias{anyRequired}
\alias{gedRequired}
\alias{isRequired}
\alias{onlyRequired}
\title{Utility Functions for GED Algorithms}
\usage{
  gedRequired(x)

  isRequired(type, x)

  onlyRequired(type, x)

  anyRequired(x)
}
\arguments{
  \item{x}{an object that defines a GED algorithm, i.e. a
  \code{GEDalgorithm} object, an access key, a registry
  entry, etc..}

  \item{type}{a character string giving the type of input
  to test, i.e. either \code{'Basis'}, \code{'Coef'} or
  \code{'Marker'}.}
}
\description{
  \code{gedRequired} tells which input data is required by
  a given, or all GED algorithm(s).

  \code{isRequired}, \code{onlyRequired} and
  \code{anyRequired} return a single logical if \code{x} is
  not missing, or a logical vector with one element per
  registered algorithm otherwise.

  \code{onlyRequired} tells if a given input data is the
  only required data for a GED algorithm to run.

  \code{anyRequired} tells if any input data is required at
  all to run a GED algorithm.
}
\details{
  If \code{x} is not missing \code{gedRequired} returns a
  logical vector with one element per type of possible
  input: \code{'Basis'}, \code{'Coef'}, \code{'Marker'}. If
  \code{x} is missing, then it returns a matrix of these
  vectors, with one row per algorithm registered in the
  internal registry.
}
\examples{
# which algorithm requires what
gedRequired()
# ask for a given algorithm
gedRequired('qprog')
}

