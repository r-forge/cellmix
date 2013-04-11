\name{gedAlgorithm}
\alias{gedAlgorithm}
\alias{gedAlgorithmInfo}
\title{Managing Gene Expression Deconvolution Algorithms}
\usage{
  gedAlgorithm(key, error = TRUE, exact = FALSE,
    all = FALSE, ...)

  gedAlgorithmInfo(show = TRUE, all = !show)
}
\arguments{
  \item{key}{algorithm access key, as a single character
  string. If missing the function returns the list of
  registered keys, as a character vector.}

  \item{error}{a logical that indicates whether an error
  should be thrown if the key is not found in the registry
  or match multiple lists. If \code{FALSE} then function
  returns \code{NULL} if the key is not found or \code{NA}
  in case of multiple matches.}

  \item{exact}{logical that indicates if one should use
  exact matching or partial matching to match the provided
  access key against all registered keys.}

  \item{all}{logical that is only used when \code{key} is
  missing and indicates if all registered keys should be
  returned, including the internal one -- whose key starts
  with a '.'.}

  \item{...}{extra arguments used internally, not to be
  used by the end user.}

  \item{show}{logical that indicates if the registry object
  should be printed (\code{FALSE}) or only returned
  invisibly (\code{FALSE}).}
}
\description{
  \code{gedAlgorithm} provides access to the registered
  gene expression deconvolution algorithms.

  \code{gedAlgorithmInfo} prints information about the
  registered gene expression deconvolution or returns --
  invisibly -- the complete algorithm registry, as a
  \code{registry} object.
}
\examples{
gedAlgorithm()
gedAlgorithm('ssKL')
gedAlgorithm('csSAM')
# show algorithms and properties
gedAlgorithmInfo()
class(gedAlgorithmInfo(show=FALSE))
}

