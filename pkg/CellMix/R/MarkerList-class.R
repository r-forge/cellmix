# Class for Marker Lists
# 
# Author: Renaud Gaujoux
# Created: 26 Jul 2012
###############################################################################

#' @include AllGenerics.R 
NULL

#' Class for Marker Gene Lists
#' 
#' This is the main class used to store sets of marker genes, for use in gene 
#' expression deconvolution methods.
#' It serves as a light version of class \code{\linkS4class{MarkerSetCollection}}, 
#' assuming that all marker identifiers are of the same type (e.g., entrez gene, unigene 
#' or probe ids), and are relative to the same platform, organism, etc.. (see details).
#' 
#' The objective of the \code{MarkerList} class is to simplify the structure and 
#' processing of gene lists, compared to handling \code{MarkerSetCollection} objects.
#' Conversion methods between these two classes are provided.  
#' 
#' The class is essentially a named list, in which each element contains data 
#' about a set of marker features/genes, e.g., from a given cell-type.
#' Marker features are assumed to be exclusive to each set, i.e. they appear in only one 
#' element of the list.
#' In addition to their identifiers, markers can be associated with numeric values, 
#' e.g., corresponding to p-values or specificity scores, or integer, which are then 
#' interpreted as indexes relative to some expression data.
#' 
#' It contains the same slots as \code{\linkS4class{GeneSet}}, which are used when 
#' converting \code{MarkerList} objects into \code{MarkerSetCollection} objects, 
#' to fill the slots of all gene set.
#' For a description of each slot, please see the documentation for 
#' \code{\linkS4class{GeneSet}}.
#' 
#' @section Heterogeneous marker lists:
#' Due to the structure of \code{MarkerList} objects, all gene identifiers must share
#' the same set of characteristics.
#' However, it is possible to associate a \code{MarkerList} object with multiple annotation
#' packages, which will be correctly dealt with if all the list's identifiers are found in only 
#' one of them.
#' For example the registered maker list \sQuote{Abbas} is associated with annotation packages
#' \code{hgu133a.db} and \code{hgu133b.db}. 
#' 
#' @import GSEABase
#' @author Renaud Gaujoux
#' @rdname markers
#' @export
setClass('MarkerList', contains='namedList',
		representation = representation(
			geneIdType = 'GeneIdentifierType',
			setName = 'ScalarCharacter',
			setIdentifier = 'ScalarCharacter',
			shortDescription = 'ScalarCharacter',
			longDescription = 'ScalarCharacter',
			organism = 'ScalarCharacter',
			pubMedIds = 'character',
			urls = 'character',
			contributor = 'character',
			version = 'Versions',
			creationDate = 'character',
			collectionType = 'CollectionType'
		), 
		prototype = prototype(
			setName = new("ScalarCharacter", NA),
			setIdentifier = new("ScalarCharacter", NA),
			geneIdType = new("NullIdentifier"),
			version = new("Versions", "0.0.1"),
			collectionType = new("NullCollection")
		),
		validity = function(object){
			
			# nothing else to be checked for empty lists
			if( length(object) == 0L ) return(TRUE)
			
			.pbdesc <- function(i, ...){
				if( !is.null(names(object)) ) str_out(names(object)[i], ...) 
				else str_out(i, ..., quote=FALSE)
			}
			# check type of elements: must all be character, integer, numeric or logical
			ok <- sapply(object, function(x) is.vector(x) && !is.list(x))
			if( !all(ok) ){
				ipb <- which(!ok)
				if( length(ipb) > 0L )
					return(str_c("Wrong type for element(s) [",.pbdesc(ipb),"]: all elements must be vectors [", class(object[[1L]]), "]."))
			}
			
			# check type of elements: must all be of the same type
			cl.ok <- sapply(object, function(x) identical(class(x), class(object[[1L]])))
			if( !all(cl.ok) ){
				ipb <- which(!cl.ok)
				if( length(ipb) > 0L )
					return(str_c("Wrong type for element(s) [",.pbdesc(ipb),"]: all elements must be vectors of the same type [", class(object[[1L]]), "]."))
			}
			
			# check presence of names for non integer lists
			has.names <- sapply(object, function(x) !is.null(names(x)))
			lelmt <- sapply(object, length)
			if( hasValues(object) ){
				ipb <- which(!has.names & lelmt)
				if( length(ipb) > 0L ){
					return(str_c("Element(s) [",.pbdesc(ipb),"] have no names: lists of non integer elements must have named elements."))
				}
			}else if( !all(has.names || !lelmt) && !all(!has.names) ) # all elements or none should be named
				return("Invalid list of non-valued elements: all elements or none of them must have names.")
			
			# check for duplicated in names
			dup <- sapply(object, function(x){
						res <- 
								if( hasValues(x) || is.logical(x) ) names(x)[anyDuplicated(names(x))]
								else if( is.character(x) ) x[anyDuplicated(x)]
								else character()
						if( length(res) == 0L ) NA else res
					})
			if( length(ipb <- which(!is.na(dup))) > 0L ){
				return(str_c(length(ipb), " top element(s) [", .pbdesc(ipb), "] have duplicated names "
								, "[", str_out(dup[ipb]) ,"]: all named elements must have unique names."))
			}
			
			# all is good
			TRUE
		}
)

#' Factory Method for Marker Lists
#' 
#' \code{MarkerList} is an S4 generic function that provides a convenient interface 
#' to create \code{\linkS4class{MarkerList}} objects from a variety of input formats,
#' such as plain named lists, factors, matrices or \code{\linkS4class{ExpressionSet}} 
#' objects, text files, etc..
#' 
#' Except for its method \code{MarkerList,matrix}, this function does not try to 
#' infer markers from the input data, but only reorganise groups of markers that 
#' are already clearly defined in the input data.
#' 
#' @param object input object from which marker data are extracted.
#' @param ... extra arguments to allow extension.
#' See each method's description for more details.
#' 
#' @inline
#' @export 
setGeneric('MarkerList', function(object, ...) standardGeneric('MarkerList'))

#' Accessing Data in Marker Lists
#' 
#' The S4 generics \code{geneValues} and \code{geneValues<-} respectively return 
#' and set the values associated with each gene in a gene set. 
#' 
#' @param object an object from/which to get/set data
#' @param value replacement value
#' @param ... extra arguments to allow extension. 
#' See each method's description for more details. 
#' 
#' @inline
#' @export
#' @rdname MarkerList-data
setGeneric("geneValues", function(object, ...) standardGeneric('geneValues'))
#' @inline
#' @rdname MarkerList-data
#' @export
setGeneric("geneValues<-", function(object, ..., value) standardGeneric('geneValues<-'))

#' \code{geneIds} returns all gene identifiers in a \code{MarkerList} 
#' object as a standard list.
#'
#' @rdname MarkerList-data
#' @inline
#' @export
setGeneric('geneIds', package='GSEABase')
#' @rdname MarkerList-data
#' @inline
#' @export
setGeneric('geneIds<-', package='GSEABase')

#' @rdname MarkerList-data
#' @inline
#' @export
setGeneric('geneIdType', package='GSEABase')
#' @rdname MarkerList-data
#' @inline
#' @export
setGeneric('geneIdType<-', package='GSEABase')

#' \code{annotation} and \code{annotation<-} respectively extract and set 
#' the annotation string from a \code{\link{MarkerList}} object, 
#' i.e. the name of Bioconductor annotation package(s) relevant to convert marker identifiers 
#' into other IDs.
#' 
#' @export
#' @rdname MarkerList-data
#' @inline
setGeneric('annotation', package='Biobase')
#' @export
#' @rdname MarkerList-data
#' @inline
setGeneric('annotation<-', package='Biobase')


#' @rdname MarkerList-data
#' @importMethodsFrom GSEABase details
#' @export
#' @inline
setGeneric('details', package='GSEABase')

#' @importMethodsFrom GSEABase mapIdentifiers
setGeneric('mapIdentifiers', package='GSEABase')

#' Identifying Gene or Probe ID Type
#' 
#' The S4 generic \code{idtype} automatically determine the type 
#' of gene/feature identifiers stored in objects, based on a combination of 
#' regular expression patterns and test functions.
#' 
#' It uses a heuristic based on a set of regular expressions and functions 
#' that uniquely match most common types of identifiers, such as Unigene, 
#' entrez gene, Affymetrix probe ids, Illumina probe ids, etc..
#' 
#' @param object an R object that contains the gene identifiers whose type is to 
#' be determined.
#' @param ... extra argument to allow extension, generally passed down to 
#' \code{idtype,character-method}.
#' See each method's description for more details.
#' 
#' @export
#' @inline
setGeneric('idtype', function(object, ...) standardGeneric('idtype'))


#' @rdname MarkerList-data
#' @export
#' @inline
setGeneric('nmf', package='NMF')
