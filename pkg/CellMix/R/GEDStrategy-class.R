# Class for GED algorithms
# 
# Author: Renaud Gaujoux
# Creation: 21 Jan 2012
###############################################################################

#' @include registry-algorithms.R
NULL

#' Class for GED Algorithms
#' 
#' 
setClassRegistry(GEDalgorithm_registry, 'GEDStrategy'
				, prototype = prototype(reqBasis=FALSE, reqCoef=FALSE, reqMarker=FALSE, cite=''))

#' @rdname GEDStrategy-class
setMethod('show', 'GEDStrategy', 
	function(object){
		cat("<Object of class: ", class(object), ">\n", sep='')
		cat("Key: ", object@key, "\n", sep='')
		req <- gedRequired(object)
		req <- names(req)[req]
		cat("Require: ", str_out(sub("^req", "", req), NA, quote=FALSE), "\n", sep='')
		if( nchar(object@cite) > 0L ){
			ref <- object@cite
			if( !grepl(" ", ref) )
				ref <- packageReference(object@cite, short=TRUE)
			cat("Reference: ", ref, "\n", sep='')
		}
	}
)

setGeneric('name', package='NMF')
setMethod('name', 'GEDStrategy', function(object) object@key )

setOldClass('GEDalgorithm_entry')
#' Factory Methods for GEDStrategy Objects
#' 
#' @param object object from which of a \code{GEDalgorithm} is created
#' @param ... extra arguments to allow extension  
#'
#' @inline 
#' @export
setGeneric('GEDStrategy', function(object, ...) standardGeneric('GEDStrategy'))
#' This method serves as a standard constructor and is equivalent to 
#' \code{new('GEDalgorithm', ...)}, meaning that all arguments must be named and 
#' correspond to slots in \code{\linkS4class{GEDStrategy}}.
#' 
setMethod('GEDStrategy', 'missing', function(object, ...) new('GEDStrategy', ...))
#' Method to create a \code{GEDalgorithm} object from its entry in the internal registry.
setMethod('GEDStrategy', 'GEDalgorithm_entry', function(object) do.call('GEDStrategy', object))
#' Constructor-Copy method, which creates an \code{GEDalgorithm} object from a template
#' and ensure that the object has a key that does not conflict with registered algorithms.
#' 
#' @param key if not missing, it must be a character string, which is 
#' used as unique key for the new GED algorithm object.
#' In particular, it cannot be a key already used by one of the registered 
#' algorithms, except if \code{object} and \code{key} are the only arguments
#' and if the key matches the key of \code{object} itself, which simply returns 
#' the object unchanged -- and is somehow not a very useful call.
#' If missing, a unique key is generated using a digest hash of the result from 
#' \code{\link{tempfile}()}.
#' 
setMethod('GEDStrategy', 'GEDStrategy', 
	function(object, key, ...){
		
		if( nargs() == 1L ) return(object)
		
		if( missing(key) ){ # generate a key if necessary
			key <- str_c('GEDStrategy_', digest(tempfile()))
		}else if( name(object) == key ){
			if( nargs() > 2L )
				stop('GEDStrategy - Cannot create a GEDStrategy object using the same key as its template [', key, "']")
			return(object)
		}else{
			meth <- gedAlgorithm(key, error=FALSE, exact=TRUE)
			if( is(meth, 'GEDStrategy') )
				stop('GEDStrategy - Cannot create GEDStrategy with key ', key, "': another algorithm already exists under the same key.")
		}
		
		# create an object based on the template
		new('GEDStrategy', meth, key=key, ...)
	}
)
#' This method creates a \code{GEDalgorithm} object from an access key, 
#' which allows to make objects using any registered algorithm as a 
#' template.
#' 
#' The access key must be valid, i.e. that it must correspond exactly 
#' to a registered algorithm.
#' If other arguments are passed, these must correspond to slots in 
#' \code{\linkS4class{GEDStrategy}}, and will overwrite the values 
#' defined by the registered algorithm. 
#' 
setMethod('GEDStrategy', 'character', 
	function(object, ...){
		
		# check if an algorithm already exists with the key
#		meth <- gedAlgorithm(object, error=nargs() == 1L, exact=TRUE)	
		meth <- gedAlgorithm(object, error=TRUE, exact=TRUE)
		
		# create/update a strategy object using provided slots if necessary 	
		if( !is(meth, 'GEDStrategy') ) new('GEDStrategy', key=object, ...)
		else if( nargs() == 1L ) meth
		else GEDStrategy(meth, ...)
		
	}
)

#' This method creates a \code{GEDalgorithm} object from a plain 
#' R function.
#' 
#' @param name descriptive name of the algorithm. 
#' A temporary name is generated if missing.
#' 
setMethod('GEDStrategy', 'function', 
	function(object, ..., name=NULL){
		
		if( is.null(name) ) name <- basename(tempfile('gedAlgorithm_'))
		GEDStrategy(key=name, algorithm=object)
		
	}
)

#' Utility Functions for GED Algorithms
#'  
#' \code{gedRequired} tells which input data is required by a given, 
#' or all GED algorithm(s).
#' 
#' If \code{x} is not missing \code{gedRequired} returns a logical vector with one 
#' element per type of possible input:
#' \code{'Basis'}, \code{'Coef'}, \code{'Marker'}.
#' If \code{x} is missing, then it returns a matrix of these vectors, with one 
#' row per algorithm registered in the internal registry.
#' 
#' @param x an object that defines a GED algorithm, i.e. a
#' \code{GEDalgorithm} object, an access key, a registry entry, etc.. 
#' 
#' @rdname GEDalgorithm-utils
#' @export
#' @examples 
#' 
#' # which algorithm requires what
#' gedRequired()
#' # ask for a given algorithm
#' gedRequired('qprog')
#' 
gedRequired <- function(x){
	
	if( missing(x) ) x <- gedAlgorithm()
	if( is.character(x) && length(x) > 1L ) return(t(sapply(x, gedRequired)))
	
	req <- 
	if( is(x, 'GEDalgorithm_entry') ){
		req <- grep("^req", names(x), value=TRUE)
		unlist(x[req])
	}else{
		# load the algorithm
		x <- GEDStrategy(x)
		req <- grep("^req", slotNames(x), value=TRUE)
		vapply(req, function(s) slot(x, s), logical(1))		
	}

	setNames(req, sub("^req","",names(req)))
}

#' \code{isRequired} tells if a given input data type is required by a GED algorithm.
#' 
#' \code{isRequired}, \code{onlyRequired} and \code{anyRequired} return a single 
#' logical if \code{x} is not missing, or a logical vector with one element per 
#' registered algorithm otherwise. 
#' 
#' @param type a character string giving the type of input to test,
#' i.e. either \code{'Basis'}, \code{'Coef'} or \code{'Marker'}.  
#' 
#' @export
#' @rdname GEDalgorithm-utils
isRequired <- function(type, x){
	
	if( missing(x) ){
		req <- gedRequired()
		req[,type]
	}else{
		req <- gedRequired(x)
		req[type]
	}
}

#' \code{onlyRequired} tells if a given input data is the only required data for 
#' a GED algorithm to run. 
#' 
#' @export
#' @rdname GEDalgorithm-utils
onlyRequired <- function(type, x){
	
	if( missing(x) ){
		req <- gedRequired()
		itype <- which(colnames(req)==type)
		ro <- !apply(req[,-itype,drop=FALSE], 1L, any)
		setNames(ro, rownames(req))
	}else{
		req <- gedRequired(x)
		itype <- which(names(req)==type)
		req[itype] && !any(req[-itype])
	}
	
}

#' \code{anyRequired} tells if any input data is required at all to run a GED 
#' algorithm.
#' 
#' @export
#' @rdname GEDalgorithm-utils
anyRequired <- function(x){
	if( missing(x) ){
		req <- gedRequired()
		ro <- apply(req, 1L, any)
		setNames(ro, rownames(req))
	}else{
		any(gedRequired(x))
	}
}
