\name{npure}
\alias{lpure}
\alias{mixedSamples}
\alias{npure}
\alias{pureSamples}
\alias{wpure}
\title{Acessing Pure or Mixed Sample Data}
\usage{
  npure(object)

  lpure(object)

  wpure(object)

  pureSamples(object, drop = TRUE)

  mixedSamples(object, drop = TRUE)
}
\arguments{
  \item{object}{an \code{\linkS4class{ExpressionMix}}
  object.}

  \item{drop}{logical that indicates if empty levels should
  be removed from phenotypic factor covariates of the
  returned \code{ExpressionMix} object.}
}
\description{
  \code{npure} returns the number of pure samples.

  \code{lpure} returns the pure samples as a logical
  vector.

  \code{wpure} returns the indexes of the pure samples.

  \code{pureSamples} returns the data for the pure samples
  only.

  \code{mixedSamples} returns the data for the mixed
  samples only.
}

