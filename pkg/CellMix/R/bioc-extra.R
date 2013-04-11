# Bioconductor extensions/utility functions
# 
# Author: Renaud Gaujoux
# Creation: 10 Jan 2012
###############################################################################

#' @include AllGenerics.R 
#' @include MarkerList-class.R
NULL

#' Numeric Computations on ExpressionSet objects 
#' 
#' @description
#' The \pkg{CellMix} defines some generics and methods to apply numeric transformations 
#' to \code{ExpressionSet} objects, which is convenient when working on gene expression 
#' deconvolution algorithms, where scale (log/linear) may matter. 
#'  
#' \code{log} log-transforms the expression matrix of \code{\link{ExpressionSet}} objects.
#' 
#' @param x an \code{ExpressionSet} object.
#' @param ... extra arguments passed to subsequent calls, usually of the corresponding  
#' method in the \pkg{stats} package.
#' 
#' @rdname ExpressionSet-compute
#' @export
setMethod('log', 'ExpressionSet'
		, function(x, ...){
			exprs(x) <- log_transform(exprs(x), ...)
			x
		}
)

#' \code{expb} applies an entry-wise exponential transformation to the expression matrix 
#' of \code{\link{ExpressionSet}}, in a specific base.
#' 
#' @export
#' @rdname ExpressionSet-compute
#' @inline
setGeneric('expb', function(x, ...) standardGeneric('expb') )
#' @param base log base to use.
#' @export
setMethod('expb', 'matrix'
		, function(x, base=exp(1)){	
			if( !missing(base) ) exp(exprs(x) * log(base)) else exp(exprs(x))
		}
)
#' @export
setMethod('expb', 'ExpressionSet'
		, function(x, ...){	
			exprs(x) <- expb(exprs(x), ...) 
			x
		}
)
#' \code{exp} is equivalent to \code{expb(x, 1)}. 
#' 
#' @export
#' @rdname ExpressionSet-compute
setMethod('exp', 'ExpressionSet', function(x) expb(x) )

#' \code{range} computes the range of expression values from 
#' \code{\link{ExpressionSet}} objects.
#' 
#' @export
#' @rdname ExpressionSet-compute
setMethod('range', 'ExpressionSet'
		, function(x, ..., na.rm = FALSE){
			range(exprs(x), ..., na.rm=na.rm)
		}
)

#' \code{quantile} computes the range of expression values in 
#' \code{\link{ExpressionSet}} objects.
#' 
#' @S3method quantile ExpressionSet
#' @rdname ExpressionSet-compute
quantile.ExpressionSet <- function(x, ...){
	quantile(exprs(x), ...)
}

#' Combining Expression Matrices
#' 
#' The method \code{cbind.ExpressionSet} combines sample expression 
#' profiles from multiple \code{\link{ExpressionSet}} or matrix objects.
#' 
#' The expression matrices must be exactly of the same dimensions.
#' For the end result to be meaningful, one probably wants the row names 
#' to match as well, i.e. both object contain the same features, in the same order. 
#' However no check is done for this.
#' 
#' Note that the returned \code{ExpressionSet} object has no sample or feature 
#' annotations.
#' 
#' @param ... series of \code{ExpressionSet} and/or matrix objects.
#' @inheritParams base::cbind
#' 
#' @return an \code{ExpressionSet} object
#' 
#' @S3method cbind ExpressionSet
#' @export
cbind.ExpressionSet <- function(..., deparse.level = 1){
	objects <- list(...)
	nobjects <- length(objects)
	out <- exprs(objects[[1]])
	#other <- names(objects[[1]]$other)
	if (nobjects > 1){ 
		lapply(2:nobjects, function(i){
					# expression values
					o <- objects[[i]]
					out <<- cbind(out, if( is(o, 'ExpressionSet') ) exprs(o) else o )
					NULL
				})
	}
	ExpressionSet(out)	
}


#' Mapping Bioconductor Annotation packages to GEO GPL Identifiers
#' 
#' \code{gpl2bioc} and \code{bioc2gpl} maps a vector of GEO GPL ids into 
#' Bioconductor annotation packages, and vice versa. 
#' 
#' @param x character vector of GPL ids or Bioconductor packages 
#' to map.
#' 
#' @return a character vector
#' 
#' @rdname Bioc2GPL
#' @export
#' @examples
#' gpl2bioc() 
#' gpl2bioc("GPL96")
#' gpl2bioc(c("GPL96", "GPL570", '123456'))
#' gpl2bioc("96")
#' 
gpl2bioc <- function(x){
	GPL2bioc <- ldata('GPL2bioc')
	map <- setNames(biocann_pkgname(GPL2bioc$bioc_package), GPL2bioc$gpl)
	map <- map[map != '']
	if( missing(x) ) return(map)
	
	res <- setNames(rep(NA, length(x)), x)
	rmgpl <- function(x) sub("^GPL", '', x)
	i <- match(rmgpl(x), rmgpl(names(map)))
	ok <- !is.na(i)
	res[ok] <- map[i[ok]]
	res
}
#' @rdname Bioc2GPL
#' @export
#' @examples 
#' bioc2gpl()
#' bioc2gpl('hgu133a')
#' bioc2gpl(c('hgu133a', 'hgu133plus2', 'abcd'))
bioc2gpl <- function(x){
	
	# get map
	map <- gpl2bioc()
	map <- revmap(map[map != ''])
	if( missing(x) ) return( map )
	
	res <- setNames(rep(NA, length(x)), x)
	i <- match(biocann_pkgname(x), biocann_pkgname(names(map)))
	ok <- !is.na(i)
	res[ok] <- map[i[ok]]
	res
}

#' Map Between Bioconductor Annotation Packages and GEO GPL Identifiers
#' 
#' The data \code{GPL2bioc} contains a data.frame that was built using the function 
#' \code{getBiocPlatformMap} from the \pkg{GEOmetadb} package.
#' 
#' @name GPL2bioc
#' @docType data
NULL


