# Definition of class MarkerSetCollection
# 
# Author: Renaud Gaujoux
# Creation: 13 Jul 2012
###############################################################################

#' @include AllGenerics.R 
NULL

library(GSEABase)

#' Classes for Marker Gene Set Collections
#' 
#' This class is an extension of the \code{\linkS4class{GeneSetCollection}} that 
#' ensures that gene/feature sets are disjoint, i.e. that each gene/feature 
#' appear in only one set. e.g., it is used to store cell-specific marker genes
#' that are used in gene expression deconvolution methods (cf. \code{\link{ged}}).
#' 
#' The main difference with \code{\linkS4class{GeneSetCollection}} is that  
#' the validity method additionally ensures that there is no duplicated gene 
#' across all sets in the collection.
#'
#' @export  
setClass('MarkerSetCollection', contains='GeneSetCollection'
	, validity= function(object){
		
		msg <- NULL
		# get all ids
		ids <- unlist(sapply(object, geneIds))
		
		# check for duplicates
		# NB: we know that there is no duplicates within each set since this is
		# checked in validity,GeneSet
		if( length(w <- which(duplicated(ids))) ){
			msg <- c(msg, str_c("MarkerSetCollection objects can only contain disjoint gene sets:"
						, " found ", length(w), " duplicated ids :", str_out(ids[w], 10)))
		}
		
		# return error if any and TRUE otherwise
		if( !is.null(msg) ) paste(msg, collapse="\n  + ")
		else TRUE
	}
)

#' Factory Methods for Marker Gene Set Collections
#'
#' @param object source object from which the object must be created
#' @param ... extra arguments to allow extension and passed to subsequent 
#' calls.
#' See details in each method's description.
#'  
#' @export 
setGeneric('MarkerSetCollection', function(object, ...) standardGeneric('MarkerSetCollection'))

