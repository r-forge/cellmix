\name{cellmix.options}
\alias{cellmix.getOption}
\alias{cellmix.options}
\alias{cellmix.printOptions}
\alias{cellmix.resetOptions}
\title{Package Specific Options for CellMix}
\usage{
  cellmix.options(...)

  cellmix.getOption(x, default = NULL)

  cellmix.resetOptions(..., ALL = FALSE)

  cellmix.printOptions()
}
\arguments{
  \item{...}{option specifications. For
  \code{cellmix.options} this can be named arguments or a
  single unnamed argument that is a named list (see
  \code{\link{options}}.

  For \code{cellmix.resetOptions}, this must be the names
  of the options to reset.}

  \item{ALL}{logical that indicates if options that are not
  part of the default set of options should be removed.}

  \item{x}{a character string holding an option name.}

  \item{default}{if the specified option is not set in the
  options list, this value is returned.  This facilitates
  retrieving an option and checking whether it is set and
  setting it separately if not.}
}
\description{
  Functions to get/set options that are specific to the
  \pkg{CellMix} package.

  \code{cellmix.options} sets/get single or multiple
  options, It behaves in the same way as
  \code{\link[base]{options}}.

  \code{cellmix.getOption} returns the value of a single
  option. It behaves in the same way as
  \code{\link[base]{getOption}}.

  \code{cellmix.resetOptions} reset all options to their
  default values.

  \code{cellmix.printOptions} prints all options along with
  their default values, in a relatively compact way.
}
\details{
  All these option management functions are generated using
  the function \code{\link[pkgmaker]{setupPackageOptions}}
  from the \pkg{pkgmaker} package.
}
\examples{
# show all specific options
cellmix.printOptions()

# get some options
cellmix.getOption('verbose')
# set new values
cellmix.options(verbose=TRUE)
cellmix.printOptions()

# reset to default
cellmix.resetOptions()
cellmix.printOptions()
}

