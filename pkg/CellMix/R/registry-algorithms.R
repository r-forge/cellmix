# Registry for deconvolution algorithms
# 
# Author: Renaud Gaujoux
# Creation: 16 Jan 2012
###############################################################################

#' @include registry.R
NULL

library(NMF)

##############################################################################
# REGISTRY DEFINITION
##############################################################################
#' Registry Object for GED Algorithms
#' 
#' Registry object of class \code{\link[registry:regobj]{registry}} for storing data on 
#' gene expression deconvolution algorithms.
#' The data are stored as plain lists but are wrapped into \code{\linkS4class{GEDStrategy}} 
#' objects when retrieved with the factory function \code{\link{GEDStrategy}}.
#' 
#' @keywords internal
GEDalgorithm_registry <- registry::registry(registry_class = "GEDalgorithm_registry"
		, entry_class = "GEDalgorithm_entry")

#' @importFrom registry match_partial_ignorecase
GEDalgorithm_registry$set_field("key", type="character"
		, is_key = TRUE
		, index_FUN = match_partial_ignorecase
		, validity_FUN = checkKey)

# Description
GEDalgorithm_registry$set_field("description", type="character", default = '')
# Algorithm
GEDalgorithm_registry$set_field("algorithm", type="function", is_mandatory = TRUE)
# Data requirement
GEDalgorithm_registry$set_field("reqBasis", type='logical', is_mandatory = TRUE)
GEDalgorithm_registry$set_field("reqCoef", type='logical', is_mandatory = TRUE)
GEDalgorithm_registry$set_field("reqMarker", type='logical', is_mandatory = TRUE)
#GEDalgorithm_registry$set_field("canfit", type='function', default = function(...) TRUE)
# Iterations
GEDalgorithm_registry$set_field("maxIter", type="integer", is_mandatory=TRUE)
# Default arguments
GEDalgorithm_registry$set_field("defaults", type="list", default=list())
# Reference
GEDalgorithm_registry$set_field("cite", type="character", default='')

# register using pkgmaker registry features
GEDalgorithm_registry <- setPackageRegistry('algorithms', GEDalgorithm_registry
						, description = 'Gene expression deconvolution algorithms'
						, entrydesc = 'GED algorithm'
						, overwrite = FALSE)
##############################################################################
##############################################################################

.altCite <- function(ref, latex=FALSE, citet=FALSE, ...){
	ref <- as.character(ref)
	if( length(bib <- grep("^[^ ]+$", ref)) > 0L ){
		ref[bib] <- 
		if( !latex ) packageReference(ref[bib], ...)
		else paste("\\cite", if(citet) 't' else '', "{", ref[bib], "}", sep='') 
	}
	ref
}

#' @S3method format GEDalgorithm_entry
format.GEDalgorithm_entry <- function(x, ..., cite=FALSE){
	
	req <- gedRequired(x)
	c(Key = x$key, Description = x$description
			, req
			, Iter = x$maxIter
			, Reference = if( cite ) .altCite(x$cite, short=TRUE) else x$cite
	)
}

#' @S3method format GEDalgorithm_registry
format.GEDalgorithm_registry <- function(x, ..., all=TRUE){
	rec <- x$get_entries()
	if( !all ){ # do not show internal methods 
		rec <- rec[grepl("^[^.]", names(rec))]
	}
	data.frame(t(sapply(rec, base::format, ...))[,-1,drop=FALSE]
				, stringsAsFactors = FALSE
				, check.names=FALSE)	 
}

#' @S3method print GEDalgorithm_registry
print.GEDalgorithm_registry <- function(x, ...) print(format(x, ...))

#' @S3method xtable GEDalgorithm_registry
#' @importFrom xtable xtable
xtable.GEDalgorithm_registry <- function(x, citet=FALSE, all=FALSE, ...){
	d <- format(x, cite=FALSE, all=all)
	# citations
	d$Reference <- .altCite(d$Reference, latex=TRUE, citet=citet)
	d$Description <- str_c(d$Description, ifelse(nchar(d$Reference), paste0(" (", d$Reference, ")"), ''))
	d$Reference <- NULL
	# required data
	req <- gedRequired(rownames(d))
	req <- t(apply(req, 1L, function(v) ifelse(v, "\\checkmark", "-")))
	d[, colnames(req)] <- req
	# iterations
	d$Iter <- ifelse(d$Iter=='1', '-', d$Iter)
	xtable::xtable(d, ...)
}

setOldClass('GEDalgorithm_registry')
setMethod('keys', signature(x='GEDalgorithm_registry')
	, function(x, keytype){
		setNames(x$get_field_entries('key'), NULL)
	}
)

#' Managing Gene Expression Deconvolution Algorithms
#' 
#' \code{gedAlgorithm} provides access to the registered gene expression 
#' deconvolution algorithms.
#' 
#' @param key algorithm access key, as a single character string. 
#' If missing the function returns the list of registered keys, as a
#' character vector.
#' @param error a logical that indicates whether an error should be thrown if 
#' the key is not found in the registry or match multiple lists.
#' If \code{FALSE} then function returns \code{NULL} if the key is not found 
#' or \code{NA} in case of multiple matches.
#' @param exact logical that indicates if one should use exact matching or 
#' partial matching to match the provided access key against all registered
#' keys.
#' @param all logical that is only used when \code{key} is missing and indicates
#' if all registered keys should be returned, including the internal one -- 
#' whose key starts with a '.'. 
#' @param ... extra arguments used internally, not to be used by the end user.
#' 
#' @export
#' @examples 
#' gedAlgorithm()
#' gedAlgorithm('ssKL')
#' gedAlgorithm('csSAM')
#' 
gedAlgorithm <- function(key, error=TRUE, exact=FALSE, all=FALSE, ...){
	
	if( missing(key) ){
		regobj <- gedAlgorithmInfo(FALSE)
		k <- keys(regobj)
		if( !all ) k <- grep("^[^.]", k, value=TRUE)
		return(k)
	}else if( is.character(key) ){
		d <- pkgreg_fetch('algorithms', key=key, error=error, exact=exact, ...)
		# exit if not found
		if( !is.list(d) ) return(d)
	}else d <- key 
	
	# create a GEDStrategy object from the entry if necessary
	GEDStrategy(d)
}

#' \code{gedAlgorithmInfo} prints information about the registered gene expression deconvolution or 
#' returns -- invisibly -- the complete algorithm registry, as a \code{registry} object.
#' 
#' @param show logical that indicates if the registry object should be printed (\code{FALSE}) 
#' or only returned invisibly (\code{FALSE}).
#' 
#' @rdname gedAlgorithm
#' @export
#' 
#' @examples 
#' # show algorithms and properties
#' gedAlgorithmInfo()
#' class(gedAlgorithmInfo(show=FALSE))
#' 
gedAlgorithmInfo <- function(show=TRUE, all=!show){
	obj <- GEDalgorithm_registry
	if( show ) print(obj, all=all)
	invisible(obj)
}

hasFormals <- function(f, args, withdefault=NA){
	has <- args %in% names(formals(f))
	if( is_NA(withdefault) ) has
	else if( has ){
		has_default <- !is.symbol(formals(f)[[args]])
		(has_default && withdefault) || (!has_default && !withdefault) 
	}
	else FALSE
}

#' Register CellMix Deconvolution Methods
#' 
#' \code{setGEDMethod} register a deconvolution algorithm 
#' in the \pkg{CellMix} registry.
#' 
#' @param key accession key.
#' @param ... extra registry fields describing the method, 
#' or arguments passed to \code{\link[pkgmaker]{pkgreg_remove}}.
#' 
#' @export
#' @rdname GEDMethods-internals
setGEDMethod <- function(key, ...){
	
	# special case for algorithms passed as character strings 
	dots <- list(...)
	if( is.character(dots$algorithm) ){
		m <- nmfAlgorithm(dots$algorithm)
		if( is(m, 'NMFStrategyIterative') 
				|| ( is(m, 'NMFStrategyFunction') && hasFormals(algorithm(m), 'maxIter') ) ){
			dots$algorithm <- eval(substitute(nmfWrapper(algorithm, maxIter=val), list(algorithm=dots$algorithm, val=dots$maxIter)))
		}else{
			dots$algorithm <- eval(substitute(nmfWrapper(algorithm), list(algorithm=dots$algorithm)))
		}
		#print(dots$algorithm)
	}
	dots$key <- key
	dots$regname <- 'algorithms'
			
	method <- do.call(setPackageRegistryEntry, dots)
	# return plain algorithm (invisible)
	return( invisible(method$algorithm) )
#	regobj <- GEDalgorithm_registry
#	if( regobj$has_entry(key) ){
#		if( !.OVERWRITE )
#		stop("GED algorithms - Could not add new entry with key '", key, "': key already exists in registry.")
#		method <- do.call(regobj$modify_entry, dots)
#	}else
#		method <- do.call(regobj$set_entry, dots)
#	
#	# return updated method as a GEDStrategy object
#	GEDStrategy(method[[1]])
}

#' \code{removeGEDMethod} removes a deconvlution algorithm 
#' from the registry.
#' Note that it does not delete any NMF algorithm that might be 
#' related to it.
#' 
#' @export
#' @rdname GEDMethods-internals
removeGEDMethod <- function(key, ...){
	pkgreg_remove('algorithms', key=key, ...)
}

#' Wrapping GED Algorithms
#' 
#' This function creates a wrapper function for calling the function 
#' \code{\link{ged}} with a given GED algorithm.
#' 
#' @param key Name of the GED algorithm to be wrapped. 
#' It should be the name of a registered algorithm as returned by \code{\link{gedAlgorithm}}, 
#' or a GED algorithm object (i.e. an instance of \code{\linkS4class{GEDStrategy}}).
#'  
#' @return a function with attribute \code{'algorithm'} set to the value of
#' \code{key}
#' 
#' @seealso \code{\link{gedAlgorithm}}, \code{\link{ged}}
#' @keywords internal
#' @export
#' 
#' @examples 
#' 
#' # wrap 'qprog' algorithm from Gong et al. (2011)
#' qp <- nmfWrapper('qprog')
#' qp
#' 
#' # test on random data
#' x <- rmatrix(100,20)
#' sig <- rmatrix(100,3)
#' res <- ged(x, sig, 'qprog', seed=12345)
#' res2 <- qp(x, sig, seed=12345)
#' nmf.equal(res, res2)
#' 
#' \dontshow{ stopifnot(nmf.equal(res, res2)) }
#' 
gedWrapper <- function(key){
	.key <- key
	f <- function(...) ged(method=.key, ...)
	attr(f, "algorithm") <- .key
	f
}

