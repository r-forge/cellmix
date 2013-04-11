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
