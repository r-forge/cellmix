\name{cellMarkers}
\alias{cellMarkers}
\alias{cellMarkersInfo}
\title{Loading Marker Lists from Registry}
\usage{
  cellMarkers(key, error = TRUE, verbose = FALSE,
    all = FALSE, ...)

  cellMarkersInfo(show = TRUE)
}
\arguments{
  \item{verbose}{logical that toggles verbose messages.}

  \item{...}{arguments passed to the marker list's loading
  function. In particular argument \code{reload=TRUE} can
  be used to clear the cache and force
  reloading/recomputing the marker list from its primary
  data file.}

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
}
\value{
  a \code{GEDmarkers_entry} object or \code{NULL} (see
  argument \code{error})
}
\description{
  Loading Marker Lists from Registry

  \code{cellMarkersInfo} prints information about the
  registered marker gene lists or return the complete
  marker list registry, as a \code{registry} object.
}
\examples{
# show marker lists and properties
cellMarkersInfo()
class(cellMarkersInfo(FALSE))
}

