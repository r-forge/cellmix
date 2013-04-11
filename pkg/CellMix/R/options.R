# Package specific options
# 
# Author: Renaud Gaujoux
# Created: 10 Dec 2012
###############################################################################

.OPTIONS <- setupPackageOptions(NAME='cellmix', RESET=TRUE,
		# default verbosity level
		verbose = FALSE
)

.IOPTIONS <- setupPackageOptions(NAME='.cellmix', RESET=TRUE,
		# wrapping ged result?
		wrap = TRUE
)

.cellmix.options <- .IOPTIONS$options
.cellmix.getOption <- .IOPTIONS$getOption

#' Package Specific Options for CellMix
#' 
#' @description
#' Functions to get/set options that are specific to the \pkg{CellMix} package.
#' 
#' \code{cellmix.options} sets/get single or multiple options,  
#' It behaves in the same way as \code{\link[base]{options}}.
#' 
#' @details
#' All these option management functions are generated using the function 
#' \code{\link[pkgmaker]{setupPackageOptions}} from the \pkg{pkgmaker} package.
#' 
#' @inheritParams base::options
#' @param ... option specifications. 
#' For \code{cellmix.options} this can be named arguments or 
#' a single unnamed argument that is a named list (see \code{\link{options}}.
#' 
#' For \code{cellmix.resetOptions}, this must be the names of the options to reset.
#' 
#' @export
#' @rdname options
#' @examples
#' 
#' # show all specific options
#' cellmix.printOptions()
#' 
#' # get some options
#' cellmix.getOption('verbose')
#' # set new values
#' cellmix.options(verbose=TRUE)
#' cellmix.printOptions()
#' 
#' # reset to default
#' cellmix.resetOptions()
#' cellmix.printOptions()
#' 
cellmix.options <- .OPTIONS$options

#' \code{cellmix.getOption} returns the value of a single option.
#' It behaves in the same way as \code{\link[base]{getOption}}.
#' 
#' @inheritParams base::getOption
#' 
#' @export
#' @rdname options
cellmix.getOption <- .OPTIONS$getOption

#' \code{cellmix.resetOptions} reset all options to their default values.
#' 
#' @param ALL logical that indicates if options that are not part of the default set 
#' of options should be removed.
#' 
#' @export
#' @rdname options
cellmix.resetOptions <- .OPTIONS$resetOptions

#' \code{cellmix.printOptions} prints all options along with their default values, 
#' in a relatively compact way.
#' @export
#' @rdname options
cellmix.printOptions <- .OPTIONS$printOptions

