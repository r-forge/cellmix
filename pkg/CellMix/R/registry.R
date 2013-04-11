# Setup of a Registry for Deconvolution methods
#
# The registry uses the `registry` package as its backbone
# http://cran.r-project/web/packages/registry/index.html
# 
# Author: Renaud Gaujoux
# Creation: 10 Nov 2011
###############################################################################

library(registry)

setGeneric('keys', package='AnnotationDbi')

#' Validity Functions of Registry Fields
#' 
#' @param x value of the field to be validated
#' @name regcheck
NULL

#' Checks if \code{x} is a valid S4 object by calling \code{\link{validObject}} 
#' on it.
#' @rdname regcheck
checkS4 <- function(x) validObject(x)
 
checkURL <- function(x){
	if( nchar(x)>0 &&  !grepl("^(http)|(ftp)://.*", x) )
		stop("registry - Invalid URL field: '", x, "'")
	TRUE
}

# access key
checkKey <- function(x){}

#' Require Existing File
#' 
#' Utility function to require a file to exist. 
#' 
#' @param file File path to check.
#' @param msg message to prepend to the error message.
#' @param sep separator between \code{msg} and general error message.
#' 
#' @keywords internal
requireFile <- function(file, msg=NULL, sep=' - '){
	if( !file.exists(file) )
		stop(if( !is.null(msg) ) str_c(msg, sep)
			, "Required file '", normalizePath(file), "' does not exist.")
	TRUE
} 


#' Advanced Fetch for Registry
#' 
#' Retrieve entries from a registry object.
#' 
#' @param regobj registry object.
#' @param key entry accession key.
#' @param msg error prefix.
#' @param error logical that indicates if an error should be thrown if the key is not found
#' or if matches multiple entries (when \var{multiple=FALSE}).
#' @param exact logical that indicates if the one looks for exact match(es).
#' @param multiple logical that indicates if the key can return multiple matches, or \code{NA} 
#' if \code{error=FALSE}, or throw an error if \code{error=TRUE}.
#' @param verbose logical that toggles verbosity.
#' 
#' @keywords internal
regfetch <- function(regobj, key, msg=NULL, error=TRUE, exact=FALSE, multiple=FALSE, verbose=FALSE){
		
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
		
	if( !is.null(msg) ) msg <- str_c(msg, ' - ')
	
	if( regobj$n_of_entries() == 0L ){
		if( error )	stop(msg, "Registry is empty: no matching entry for key ", dQuote(key), "."
						, call.=is.null(msg))
		else return(NULL)
	}
	
	d <- regobj$get_entries(key)
	
	# no entry found
	if( is.null(d) ){
		if( error ){
			stop(msg, "No matching entry for key ", dQuote(key), " in the registry."
					, "\n  ", if( !is.null(msg) ) "     "
					, "Use one of: ", str_wrap(str_out(sort(regobj$get_entry_names()), Inf), exdent=2)
					, call.=is.null(msg))
		}else return(NULL)
	}
	# multiple match
	if( length(d) > 1L ){
		# look for exact matches
		i <- which(key == vapply(d, '[[', character(1), 'key'))
		if( length(i) > 0L ) d <- d[i]
		else if( error ){
				stop(msg, "Multiple entries found for key ", dQuote(key), ": ", str_out(sort(names(d)), Inf)
					, call.=is.null(msg))
		} else if( !multiple ) return(NA)
	}
	
	# get single match
	d1 <- d[[1]]
	
	# check it is an exact match
	if( exact && d1$key != key ){
		if( error ){
			stop(msg, "No entry for '", key, "' in the registry."
				, "\n  ", if( !is.null(msg) ) "     " 
				, "Use one of: ", str_wrap(str_out(regobj$get_entry_names(), Inf), exdent=2)
				, call.=is.null(msg))
		}else return(NULL)
	}
	
	if( multiple ) d
	else d1
}
