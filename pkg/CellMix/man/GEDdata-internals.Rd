\name{GEDdata-internals}
\alias{gedDataCheck}
\alias{GEDdata-internals}
\alias{.makeData}
\alias{meanBasis}
\alias{pdataCoef}
\alias{pureCoef}
\alias{removeGEDData}
\alias{setGEDData}
\title{Internal Functions for the CellMix Dataset registry}
\usage{
  setGEDData(key, ...)

  removeGEDData(key, ...)

  meanBasis(var)

  pdataCoef(var, scale = 100)

  pureCoef(var)

  gedDataCheck(key = NULL, load = TRUE, cache = TRUE,
    save = FALSE, verbose = 3)

  .makeData()
}
\arguments{
  \item{key}{accession key. Currently, only GEO accesion
  keys are supported.}

  \item{...}{extra registry fields describing the dataset,
  or arguments passed to
  \code{\link[pkgmaker]{pkgreg_remove}}.}
}
\description{
  Internal Functions for the CellMix Dataset registry

  \code{setGEDData} register a dataset in the \pkg{CellMix}
  registry.

  \code{removeGEDData} removes a dataset from the registry.
  Note that it does not delete any local/cache file related
  to this dataset.

  \code{meanBasis} builds a basis matrix for a data entry
  from one of its phenotypic data. It tries to match the
  order of a coefficient matrix if any, or use the order of
  the variable's levels, which must have corresponding
  phenotypic variables.

  \code{pdataCoef} builds a coefficient matrix for a data
  entry from one of its phenotypic variables, whose levels
  must correspond to numeric phenotypic variables, that
  each contains proportions of one cell type that is
  present in each samples.

  \code{pureCoef} builds a coefficient matrix for a data
  entry where all samples are pure, using one of its
  phenotypic -- factor -- variable whose levels define the
  cell type of origin of each sample.

  \code{gedDataCheck} check if all or a given registered
  dataset can be loaded.

  \code{.makeData} remake some rda files shipped with
  CellMix.
}
\keyword{internal}

