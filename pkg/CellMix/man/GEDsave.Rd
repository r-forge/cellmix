\name{GEDsave}
\alias{GEDpath}
\alias{GEDsave}
\alias{GEDtmp}
\title{Internal Functions for CellMix Data Registry}
\usage{
  GEDsave(x, force = FALSE)

  GEDpath(..., create = TRUE, reset = NULL)

  GEDtmp(..., create = TRUE)
}
\arguments{
  \item{x}{GEDdata_entry object or list of GEDdata_entry
  objects. If missing then all entries are used.}

  \item{force}{logical to force saving dataset(s)}

  \item{reset}{new value for the local directory}

  \item{...}{path parts passed to \code{\link{file.path}}
  to be appended to the main path.}

  \item{create}{logical that indicates if the directory
  should be created if it does not exists.}
}
\description{
  These functions are used internally to work with the gene
  expression dataset registry within the CellMix package.

  \code{GEDpath} returns/sets the path to the local
  directory where CellMix data (e.g., cache) are stored.

  \code{GEDtmp} returns the path to the local directory
  where downloaded data from GEO are stored.
}
\details{
  \code{GEDsave} saves the primary dataset loaded in the
  registry on disk.
}
\keyword{internal}

