# Common set of S4 generics
# 
# Author: Renaud Gaujoux
# Created: 28 Nov 2012
###############################################################################

#' @include AllClasses.R
NULL

#' Extracting Feature Names
#' 
#' The \pkg{CellMix} package provides extra methods for the generic 
#' \code{\link[Biobase]{featureNames}} and \code{\link[Biobase]{sampleNames}}, 
#' that complete the original Bioconductor interface.
#' 
#' @inheritParams Biobase::featureNames
#' 
#' @rdname Bioc-generics
#' @export
#' @inline
setGeneric('featureNames', package='Biobase')

#' @rdname Bioc-generics
#' @export
setGeneric('featureNames<-', package='Biobase')

#' Returns the row names of \code{object}.
setMethod('featureNames', 'matrix', function(object) rownames(object) )

#' @rdname Bioc-generics
#' @export
#' @inline
setGeneric('sampleNames', package='Biobase')

#' @rdname Bioc-generics
#' @export
setGeneric('sampleNames<-', package='Biobase')

#' Returns the column names of \code{object}.
setMethod('sampleNames', 'matrix', function(object) colnames(object) )

#' @export
#' @rdname Bioc-generics
setGeneric('exprs', package = 'Biobase')

#' Simply returns \code{object}.
#' This method is defined so that function can seamlessly handle both matrix and 
#' \code{ExpressionSet} objects.
#' 
#' @rdname Bioc-generics
setMethod('exprs', 'matrix', function(object) object )

setGeneric('subset', package='base')

#' Enhanced Subsetting for Matrix-like Data 
#' 
#' These methods subset Matrix-like data only keeping rows 
#' whose names are in common with some other data.
#' 
#' @details 
#' The generic for these methods is taken from 
#' the \pkg{GSEABase} package because it is entirely imported
#' by \pkg{CellMix}.
#' However, a more natural definition would be to be able to use 
#' the -- currently conflicting -- version defined in the 
#' recommended \pkg{BiocGenerics} package.
#' 
#' @inheritParams base::intersect
#' @name intersect
NULL
#setGeneric('intersect', package='GSEABase')

#' This method is equivalent to \code{\link[base]{subset}(x, y)}.
#' @rdname intersect
#' @export 
setMethod('intersect', signature(x='MatrixData', y='logical'), 
	function(x, y){
		subset(x, y)
	}
)
#' Subset a matrix-like object by only keeping the rows whose 
#' names are in a given reference character vector.
#' @rdname intersect
#' @export 
setMethod('intersect', signature(x='MatrixData', y='character'), 
	function(x, y){
		intersect(x, featureNames(x) %in% y)
	}
)
#' Subset a matrix-like object by only keeping the rows whose 
#' names are in common with another matrix-like data.
#' 
#' This is a shortcut for \code{intersect(x, featureNames(y), ...)}.
#' @rdname intersect
#' @export 
setMethod('intersect', signature(x='MatrixData', y='MatrixData'), 
		function(x, y){
			intersect(x, featureNames(y))
		}
)
