\name{match_first}
\alias{match_first}
\title{Returns the index of the first match in a reference table}
\usage{
  match_first(x, table, nomatch = NA_integer_,
    incomparables = NULL)
}
\arguments{
  \item{x}{vector or \code{NULL}: the values to be
  matched.}

  \item{table}{vector or \code{NULL}: the values to be
  matched against.}

  \item{nomatch}{the value to be returned in the case when
  no match is found.  Note that it is coerced to
  \code{integer}.}

  \item{incomparables}{a vector of values that cannot be
  matched.  Any value in \code{x} matching a value in this
  vector is assigned the \code{nomatch} value.  For
  historical reasons, \code{FALSE} is equivalent to
  \code{NULL}.}
}
\description{
  This function provides a modified behaviour of
  \code{\link{match}} in the case of multiple matches,
  where it returns the index of the first match.
}
\seealso{
  \code{\link{match}}, \code{\link{charmatch}}
}

