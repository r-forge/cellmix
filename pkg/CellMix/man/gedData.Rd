\name{gedData}
\alias{gedData}
\alias{gedDataInfo}
\title{Loading Gene Expression Deconvolution Data}
\usage{
  gedData(key, error = TRUE, all = FALSE, exact = FALSE,
    with = NULL, ...)

  gedDataInfo(show = TRUE)
}
\arguments{
  \item{with}{character vector, only used when \code{key}
  is missing, that specifies some filtering criterium. It
  allows to only list datasets that contain specific data,
  which may help choosing a suitable dataset when
  developping/testing deconvolution methods: \describe{
  \item{\code{'prop'}}{cell proportions for each samples.}
  \item{\code{'sig'}}{cell-specific sigantures.}
  \item{\code{'all'}}{both proportions and cell-specific
  signatures.} \item{\code{'mixed'}}{mixed samples.}
  \item{\code{'pure'}}{pure samples.} }

  Any combination of these is allowed, e.g. \code{c('prop',
  'pure')} to list datasets that contain proportions and
  pure samples.}

  \item{show}{logical that indicates if the registry object
  should be printed (\code{FALSE}) or only returned
  invisibly (\code{FALSE}).}

  \item{key}{algorithm access key, as a single character
  string. If missing the function returns the list of
  registered keys, as a character vector.}

  \item{error}{a logical that indicates whether an error
  should be thrown if the key is not found in the registry
  or match multiple lists. If \code{FALSE} then function
  returns \code{NULL} if the key is not found or \code{NA}
  in case of multiple matches.}

  \item{all}{logical that is only used when \code{key} is
  missing and indicates if all registered keys should be
  returned, including the internal one -- whose key starts
  with a '.'.}

  \item{exact}{logical that indicates if one should use
  exact matching or partial matching to match the provided
  access key against all registered keys.}

  \item{...}{extra arguments used internally, not to be
  used by the end user.}
}
\value{
  a \code{GEDdata_entry} object or \code{NULL} (see
  argument \code{error})
}
\description{
  \code{gedData} loads entries from the CellMix dataset
  registry.

  \code{gedDataInfo} prints information about the
  registered gene expression datasets or returns --
  invisibly -- the complete dataset registry, as a
  \code{registry} object.
}
\examples{
# retrieve a dataset entry
gedData('GSE20300')

# error if the entry does not exists
try( gedData('GSE1234') )
# unless error=FALSE
gedData('GSE1234', error=FALSE)

# list datasets that contain cell proportions
gedData(with='prop')
# or pure samples
gedData(with='sig')
# or both
gedData(with=c('prop', 'sig'))
# or mixed samples
gedData(with='mixed')
# show algorithms and properties
gedDataInfo()
class(gedDataInfo(FALSE))
}

