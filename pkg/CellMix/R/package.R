# Low-level Package Functions: .onLoad, .onUnload, etc...
# 
# Author: Renaud Gaujoux
# Creation: 31 Oct 2011
###############################################################################



#' Gene Expression Deconvolution
#' 
#'
#' \tabular{ll}{
#' Package: \tab CellMix\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2011-11-01\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' @author
#' Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#'
#' Maintainer: Renaud Gaujoux \email{renaud@@cbio.uct.ac.za}
#' @name CellMix-package
#' @rdname CellMix-package
#' @docType package
#' @title Gene Expression Deconvolution
#' @keywords package
#' @bibliography ~/Documents/articles/library.bib
#' 
#' @import stringr
#' @import pkgmaker
#' @import NMF
#' 
#' @examples
#' # some technical info about the package
#' CellMix()
#' 
#' 
NULL

#' Information About the CellMix Package
#' 
#' Prints/retrieve technical information or data about the \pkg{CellMix} package.
#'  
#' @export 
#' @rdname CellMix-package
CellMix <- function(){
	
	v <- packageVersion('CellMix')
	res <- list(
			version=v
			, sversion=as.character(v)
			, userdir = GEDpath(create=FALSE)
			, registry=packageRegistry(package='CellMix')
	)
	class(res) <- 'CellMix_info'
	res
}

#' @S3method print CellMix_info
print.CellMix_info <- function(x, ...){
	cat("CellMix package info:\n")
	cat("* Version: ", x$sversion, "\n", sep='')
	cat("* Data path: ", x$userdir, "\n", sep='')
	cat("* Registries:\n")
	print(x$registry)
	cat("\n")
	invisible()
}