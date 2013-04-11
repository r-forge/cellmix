# Dataset registry definitions and functions
# 
# Author: Renaud Gaujoux
# Creation: 08 Dec 2011
###############################################################################
#' @include ExpressionMix-class.R
#' @include registry.R 
#' @include data-utils.R
NULL

#' Internal Functions for the CellMix Dataset registry
#' 
#' 
#' @name GEDdata-internals
#' @rdname GEDdata-internals
#' @keywords internal
NULL

GEDdata_registry <- registry::registry(registry_class = "GEDdata_registry"
		, entry_class = "GEDdata_entry")

# access key
GEDdata_registry$set_field("key", type="character"
		, is_key = TRUE
		, index_FUN = match_partial_ignorecase
		, validity_FUN = checkKey)
# aliases
GEDdata_registry$set_field("aliases", type="character", default = '')
# URL
GEDdata_registry$set_field("url", type="character", default = '')
# Description
GEDdata_registry$set_field("description", type="character", default = '')
# Reference
GEDdata_registry$set_field("cite", type="character", default = '')
# Annotation package
GEDdata_registry$set_field("annotation", type="character", default = '')
# Dimensions
GEDdata_registry$set_field("dim", type="integer", is_mandatory = TRUE)
# Functions to extract the signatures and mixture proportion matrices
GEDdata_registry$set_field("basis", type="function")
GEDdata_registry$set_field("coef", type="function")
# definition of the pure samples
GEDdata_registry$set_field("pure")
pureIndexes <- function(d){
	spec <- d$pure
	if( is.function(spec) ) spec
	else if( is_NA(spec) ) NA # all mixed
	else if( isTRUE(spec) ){ # all pure
		function(eset){ seq(ncol(eset)) }
	}else
		stop("Invalid definition of pure samples: must be TRUE, NA or a function [", class(spec), ']')
}
# Data itself
GEDdata_registry$set_field("env", type="list")
# Function that computes better formatted Phenotypic data 
GEDdata_registry$set_field("pdata", type="function")
# Function that computes extra Feature data 
GEDdata_registry$set_field("fdata", type="function")
# Function that filters the data 
GEDdata_registry$set_field("filter", type="function")

# register using pkgmaker registry features
GEDdata_registry <- setPackageRegistry('data', GEDdata_registry
		, description = 'Benchmarks datasets that contain cell type signatures and/or sample proportions'
		, entrydesc = 'GED dataset')
#######################################################################

esetPath <- function(object) GEDpath(paste(object$key, '.rds', sep=''))

is.GEDdata <- function(x) is(x, 'GEDdata_entry')

setOldClass('GEDdata_registry')
setMethod('keys', signature(x='GEDdata_registry')
		, function(x, keytype){
			setNames(x$get_field_entries('key'), NULL)
		}
)

# extract data source from entry key.
#' @S3method DataSource GEDdata_entry
DataSource.GEDdata_entry <- function(x, ...){
	DataSource(x$key, ...)
}

#' @S3method format GEDdata_entry
format.GEDdata_entry <- function(x, ...){
	d <- x$dim
	c(Key=x$key, Description=x$description
			, Features=d[1], Samples=d[2], Types=d[3]
			, Mixed = !isTRUE(x$pure), Pure = !is_NA(pureIndexes(x))
			, Basis=!is_NA(x$basis), Coef=!is_NA(x$coef)
			, Annotation=paste(annotation(x), collapse=', ')
			, Reference = x$cite) 	
}

#' @S3method format GEDdata_registry
format.GEDdata_registry <- function(x, ...){
	rec <- x$get_entries()
	data.frame(t(sapply(rec, base::format, ...))[,-1])	 
}

#' @S3method xtable GEDdata_registry 
#' @importFrom xtable xtable
xtable.GEDdata_registry <- function(x, citet=FALSE, ...){
	d <- format(x)
	checkmark <- function(x) ifelse(d[[x]]=="TRUE", "\\checkmark", "-")
	cn <- c('Mixed', 'Pure', 'Basis', 'Coef')
	d[cn] <- sapply(cn, checkmark)
#	d$Basis <- ifelse(d$Basis=="TRUE", "\\checkmark", "-")
#	d$Coef <- ifelse(d$Coef=="TRUE", "\\checkmark", "-")
	d$Reference <- .altCite(d$Reference, latex=TRUE, citet=citet)
	xtable::xtable(d, ...)
}

#' @S3method print GEDdata_registry
print.GEDdata_registry <- function(x, ...){	
	registry:::print.registry(x)
	print(format(x))
}


#' Loading Gene Expression Deconvolution Data
#' 
#' \code{gedData} loads entries from the CellMix dataset registry. 
#' 
#' @inheritParams gedAlgorithm
#' @param with character vector, only used when \code{key} is missing,
#' that specifies some filtering criterium.
#' It allows to only list datasets that contain specific data, which may 
#' help choosing a suitable dataset when developping/testing deconvolution 
#' methods:
#' \describe{
#' \item{\code{'prop'}}{cell proportions for each samples.}
#' \item{\code{'sig'}}{cell-specific sigantures.}
#' \item{\code{'all'}}{both proportions and cell-specific signatures.}
#' \item{\code{'mixed'}}{mixed samples.}
#' \item{\code{'pure'}}{pure samples.}
#' }
#' 
#' Any combination of these is allowed, e.g. \code{c('prop', 'pure')} to list 
#' datasets that contain proportions and pure samples.
#' 
#' @return a \code{GEDdata_entry} object or \code{NULL} (see argument \code{error}) 
#' 
#' @export
#' 
#' @examples
#' 
#' # retrieve a dataset entry
#' gedData('GSE20300')
#' 
#' # error if the entry does not exists 
#' try( gedData('GSE1234') )
#' # unless error=FALSE
#' gedData('GSE1234', error=FALSE)
#' 
#' # list datasets that contain cell proportions
#' gedData(with='prop')
#' # or pure samples
#' gedData(with='sig')
#' # or both
#' gedData(with=c('prop', 'sig'))
#' # or mixed samples
#' gedData(with='mixed')
#' 
gedData <- function(key, error = TRUE, all=FALSE, exact=FALSE, with=NULL, ...){
	
	# simple cases
	regobj <- gedDataInfo(FALSE)
	
	if( missing(key) ){
		k <- keys(regobj)
		if( !all ) k <- grep("^[^.]", k, value=TRUE)
		# filter datasets based on criteria
		if( !is.null(with) ){
			entries <- sapply(k, function(x) regobj$get_entry(x), simplify=FALSE)
			
			wset <- c('prop', 'sig', 'pure', 'mixed')
			if( any(bad <- !with %in% wset) ){
				warning("Discarding invalid filtering values ", str_out(with[bad], Inf), ".\n"
						, "  Only possible values are: ", str_out(wset, Inf))
				with <- with[!bad]
			}
			# keep data with mixture coefficient?
			if( any(c('all', 'prop', 'coef') %in% with) ){ 
				sel <- sapply(entries, function(x) !pkgmaker::is_NA(x$coef))
				entries <- entries[sel]
			}
			# keep data with mixture coefficient?
			if( any(c('all', 'sig', 'basis') %in% with) ){ 
				sel <- sapply(entries, function(x) !pkgmaker::is_NA(x$basis))
				entries <- entries[sel]
			}
			# keep data with mixed samples?
			if( any('mixed' %in% with) ){
				sel <- sapply(entries, function(x) !isTRUE(x$pure))
				entries <- entries[sel]
			}
			# keep data with pure samples?
			if( any('pure' %in% with) ){ 
				sel <- sapply(entries, function(x) !is_NA(pureIndexes(x)))
				entries <- entries[sel]
			}
			k <- names(entries)
		}
		return(k)
	}
	if( is.GEDdata(key) ) return( key )
	
	# check for an alias
	aliases <- regobj$get_field_entries('aliases', unlist=FALSE)
	if( length(wa <- which(sapply(aliases, grepl, pattern=paste0('^', key)))) ){
		if( length(wa) > 1L ){
			stop("GED dataset - Could not fetch entry: multiple aliases match query key '", key, "' [", str_out(names(aliases)[wa], Inf), "]")
		}
		# subsitute alias by key
		key <- names(aliases)[wa]
	}
	#
	
	# fetch entry in package registry
	d <- pkgreg_fetch('data', key=key, error=error, exact=exact, ...)
	
	# exit if not found
	if( !is.list(d) ) return(d)
	
	d
}

#' \code{gedDataInfo} prints information about the registered gene expression datasets or 
#' returns -- invisibly -- the complete dataset registry, as a \code{registry} object.
#' 
#' @param show logical that indicates if the registry object should be printed (\code{FALSE}) 
#' or only returned invisibly (\code{FALSE}).
#' 
#' @rdname gedData
#' @export
#' 
#' @examples 
#' # show algorithms and properties
#' gedDataInfo()
#' class(gedDataInfo(FALSE))
#' 
gedDataInfo <- function(show=TRUE){
	obj <- GEDdata_registry
	if( show ) print(obj)
	invisible(obj)
}

setOldClass('GEDdata_entry')

.onLoad.registry <- function(){
#	dir.create(GEDpath(), showWarnings = FALSE)
#	dir.create(GEDtmp(), showWarnings = FALSE)
	
#	# load NMF plugins
#	where <- packageEnv()
#	.init.nmf.plugins(where)	
}

#' Extracting Expression Data from CellMix Dataset Registry Entry.
#' 
#' Loads and extracts the expression data from a CellMix registry dataset entry, 
#' as an \code{\link{ExpressionSet}} object.
#' 
#' The first time the expression data is accessed, and if argument \code{load=TRUE}, 
#' it is loaded from the local storage if possible (cf. \code{\link{GEDsave}}), 
#' or downloaded from GEO using \code{\link{GEDdownload}}. It is then 
#' pre-processed (see below) and the resulting value (an \code{ExpressionSet} object) is 
#' cached in the GEDdata_entry environment. Subsequent calls use the cached value.
#' 
#' The data pre-processing consists in filtering and extraction/addition of 
#' phenotypic or feature annotation data. 
#' The filtering is performed by the function stored in element \code{object$filter}. 
#' The extra phenotypic data is generated by the function stored in element \code{object$pdata}. 
#' The extra feature data is generated by the function stored in element \code{object$fdata}.
#' Note that the whole pre-processing pipeline is not cached.
#' 
#' @param load logical that indicates if the primary data should be (re-)loaded, 
#' or a character string giving the path of the file from where to load the data, 
#' which will typically be the Matrix file as downloaded from GEO.
#' If \code{FALSE} the gene expression data is not loaded and the result will be 
#' \code{NULL} if the dataset has not already been loaded in a previous call.
#' 
#' If explicitly \code{TRUE} (i.e. not missing), then the data is 
#' (re-)loaded from cache, or downloaded from a public repository if it has 
#' not been previously cached.
#' In this case, the pre-processing pipeline is performed again.
#' @param filter logical that indicates if the preprocessing pipeline should
#' be applied to the "plain" data.
#' If \code{TRUE} then the processing/filtering function associated with the 
#' data (registry field \code{'filter'}) is applied.
#' A custom pre-processing function may also be passed.
#' 
#' If \code{FALSE}, no processing is performed on the data.
#' @param verbose logical/number that toggles verbose messages.
#' To better monitor how the data is generated/processed, it is recommended
#' to toggle on verbosity, and use numeric values at various levels 
#' (e.g. \code{verbose=3}) to display more detailed outputs.
#' @param save logical that indicates if the data should be saved in in R binary format 
#' in the local GED cache, for fast future loading -- possibly in a different 
#' R session.
#' @param cache logical that indicates if the data should be retrieved from the 
#' local cache if present, or re-downloaded.
#' If \code{save=TRUE} (as default), the new data will replace any previously 
#' downloaded data.
#' 
#' @return an \code{\link{ExpressionSet}} object.
#' @seealso \code{\link{ExpressionSet}}
setMethod('eset', 'GEDdata_entry',
		function(object, load = TRUE, filter = TRUE
				, verbose = FALSE, save = TRUE, cache = TRUE){
			
			# set verbosity level
			if( !missing(verbose) ){
				ol <- lverbose(verbose)
				on.exit( lverbose(ol) )
			}
			
			# check if already loaded
			env <- object$env[[1]]
			key <- object$key
			if( !identical(load, FALSE) && (is.null(env$eset) 
						|| !missing(load) 
						|| isFALSE(cache) 
						|| is.function(filter))  ){
				# load from disk
				tfile <- 'local'
				autocache <- FALSE
				p <- if( is.character(load) ) load[1] 
						else if( cache ){
							autocache <- TRUE
							tfile <- 'cache'
							esetPath(object)
						}
				if( !is.null(p) ){ # some path is specified
					if( grepl("\\.rds$", p) ) { # load binary file
						if( !autocache || file.exists(p) ){
							vmessage("Loading expression data ", key, " from ", tfile, " binary file '", p, "' ", appendLF = FALSE)
							if( !file.exists(p) ){ # error only if the file was not passed as argument
								vmessage('ERROR')
								stop("Binary data file '", p, "' does not exist.")
							}
							env$eset <- readRDS(p)
							vmessage("OK")
						} 
					} else{ # load local GEO file
						if( !file.exists(p) ){
							vmessage('ERROR')
							stop("GEO data file '", p, "' does not exist.")
						}
						vmessage("Loading expression data ", key, " from ", tfile, " GEO file '", p, "' ", appendLF = FALSE)
						env$eset <- GEOquery::getGEO(filename = p, destdir = GEDtmp())
						vmessage("OK")
					}
				}
				# download from GEO if not already loaded from disk
				if( is.null(env$eset) ){ 
					vmessage("Downloading data ", key, " from GEO '", GEDurl(object), "' ... ", appendLF = FALSE)
					env$eset <- GEDdownload(object, destdir = GEDtmp())
					vmessage("DONE")					
					if( save ){
						vmessage("Save data ", key, " in local repository ... ", appendLF=FALSE)
						GEDsave(object, force=TRUE)
						vmessage('OK')
					}
				}
				
				if( is.null(env$eset) )
					vmessage("NOTE: Expression data was not loaded")
				else if( !isFALSE(filter) ){
					
					.filter <- object$filter
					if( is.function(filter) ){
						vmessage("# NOTE: Using custom processing function")
						.filter <- filter
					}
					# filter if necessary
					if( !is_NA(.filter) ){
						if( is.list(env$eset) )
							vmessage("Processing/Filtering data [", col(sapply(env$eset, function(x) col(dim(x), sep='x')), sep=' | '),"] -> ", appendLF=FALSE)
						else
							vmessage("Processing/Filtering data [", col(dim(env$eset), sep='x'),"] -> ", appendLF=FALSE)
						fe <- .filter(env$eset)	
						# drop levels
						fData(fe) <- droplevels(fData(fe))
						if( !is(fe, 'ExpressionSet') )
							stop("Invalid object returned by `filter()` for data ", key, ": expected an ExpressionSet object [", class(fe), "]")
						if( ncol(fe) == 0 )
							stop("Invalid object returned by `filter()` for data ", key, ": no samples.")
						vmessage("[", col(dim(fe), sep='x'),"]")					
						env$eset <- fe
						rm(fe)
					}
					
					# merger function for adding extra annotations (pheno or feature)
					.merge_ann <- function(x, y, msg){
						# remove auto indentation in log
						if( verbose ) logger.options(autoindent=FALSE)
						
						# only add the variables of y that are not already in x
						i <- which(!is.element(colnames(y), colnames(x)))
						if( length(i)>0 ){
							vmessage(msg, ": ", str_wrap(col(colnames(y)[i], sep=', '), exdent=21))
							for(cov in colnames(y)[i])
								x[[cov]] <- y[[cov]]
						}
						droplevels(x)
					}
					
					# add extra phenotypic data if possible and necessary
					if( !is_NA(object$pdata) ){
						pd <- pData(env$eset) # original pdata
						pd2 <- object$pdata(env$eset) # new data
						pData(env$eset) <- .merge_ann(pd, pd2, "Adding covariate(s)")
					}
					
					# add extra feature data if possible and necessary
					if( !is_NA(object$fdata) ){
						pd <- fData(env$eset) # original pdata
						pd2 <- object$fdata(env$eset) # new data
						fData(env$eset) <- .merge_ann(pd, pd2, "Adding annotation(s)")
					}
					
					# change the annotation
					if( annotation(object)[1] != '' ) 
						annotation(env$eset) <- annotation(object) 
				}
			}
			
			env$eset
		}
)

#' Accessing Data from CellMix Dataset
#' 
#' Funtions to access data from the benchmark dataset registry entries 
#' 
#' @rdname GEDdata-access
#' @name GEDdata-access
NULL

#' \code{annotation} gets the name of the annotation package(s) relevant for the data.
#'  
#' @return a character vector, usually of length 1 but possibly longer. 
#' @export 
setMethod('annotation', 'GEDdata_entry', function(object) object$annotation)

#' \code{exprs} loads and extracts the expression matrix from a registry data entry.
#' 
#' @return a numeric matrix 
#' @rdname GEDdata-access
#' @export 
setMethod('exprs', 'GEDdata_entry',
		function(object){
			exprs(eset(object))	
		}
)

setGeneric('dims', package='Biobase')
#' \code{dims} gets the dimensions of an ExpressionMix object. It returns a 3-length integer 
#' vector, containing the number of features, samples and components respectively.
#' 
#' @rdname GEDdata-access
#' @export
setMethod('dims', 'GEDdata_entry', 
		function(object){
			object$dim
		})

setGeneric('basis', package='NMF')
#' \code{basis} access/compute the signatures and mixture Proportions.
#' 
#' @param object an object of class \code{GEDdata_entry} as returned by 
#' \code{GEDdata_registry$get_entry}.
#' 
#' @rdname GEDdata-access
#' @export
setMethod('basis', 'GEDdata_entry',
		function(object){
			f <- object$basis
			if( is_NA(f) ) NA_Matrix(0, nbasis(object))
			else if( is.function(f) ) f(object)
			else 
				stop("basis - Invalid value for field `basis` in data registry entry '", object$key, "' [", class(f), "]")
		}
)

setGeneric('coef', package='stats')
#' \code{coef} returns the mixture coefficient matrix (i.e. the proportions) available for 
#' the data.
#' 
#' @return a matrix
#' @export
#' @rdname GEDdata-access
setMethod('coef', 'GEDdata_entry',
		function(object){
			f <- object$coef
			if( is_NA(f) ) NA_Matrix(nbasis(object), 0)
			else if( is.function(f) ) f(object)
			else 
				stop("coef - Invalid value for field `coef` in data registry entry '", object$key, "' [", class(f), "]")
		}
)

setGeneric('nbasis', package='NMF')
#' \code{nbasis} returns the number of constituents signatures available in the data.
#'  
#' @return an integer
#' @export 
#' @rdname GEDdata-access
setMethod('nbasis', 'GEDdata_entry',
		function(x){
			dims(x)[3]
		}
)

#' Register CellMix Datasets
#' 
#' \code{setGEDData} register a dataset in the \pkg{CellMix} registry.
#' 
#' @param key accession key.
#' Currently, only GEO accesion keys are supported.
#' @param ... extra registry fields describing the dataset, 
#' or arguments passed to \code{\link[pkgmaker]{pkgreg_remove}}.
#' 
#' @export
#' @rdname GEDdata-internals
setGEDData <- function(key, ...){
	env <- list(new.env(parent=emptyenv()))
	res <- setPackageRegistryEntry('data', key, ..., env=env)
	invisible(res)
}

#' \code{removeGEDData} removes a dataset from the registry.
#' Note that it does not delete any local/cache file related to this 
#' dataset.
#' 
#' @export
#' @rdname GEDdata-internals
removeGEDData <- function(key, ...){
	pkgreg_remove('data', key=key, ...)
}

#' \code{meanBasis} builds a basis matrix for a data entry from one of its phenotypic data.
#' It tries to match the order of a coefficient matrix if any, or use the order of the 
#' variable's levels, which must have corresponding phenotypic variables.
#' 
#' @rdname GEDdata-internals
meanBasis <- function(var){
	function(object){
		
		# load eset
		eset <- eset(object)
		# check variable
		allv <- varLabels(eset)
		if( !var %in% allv ){
			stop("Could not build coef matrix for data entry '", object$key, "':"
					," phenotypic variable '", var, "' is not in the phenotypic data."
					, "\n  [varLabels: ", str_out(allv), ']')
		}
		# convert variable into factor if necessary
		if( !is.factor(v <- eset[[var]]) ) 
			v <- factor(v)
		
		# load coef
		co <- coef(object)
		n <- rownames(co)
		if( is.null(n) ){
			if( !is_NA(pidx <- pureIndexes(object)) ){
				v[-pidx(eset)] <- NA
				v <- droplevels(v)
			}else
				stop("Could not build mean basis for data entry '", object$key, "': no coef nor pure sample definition.")
		}else{
			# remove values not in coef
			v[!v %in% n] <- NA
			v <- droplevels(v)
			# reorder levels as in coef
			v <- factor(v, levels=n)
		}
		# compute row means in each group of samples
		rowMeansBy(exprs(eset), v)
	}
}

#' \code{pdataCoef} builds a coefficient matrix for a data entry from one of its phenotypic variables, 
#' whose levels must correspond to numeric phenotypic variables, that each contains proportions of one 
#' cell type that is present in each samples.
#' 
#' @rdname GEDdata-internals
pdataCoef <- function(var, scale=100){
	function(object){
		# load eset
		eset <- eset(object)
		# check variable
		allv <- varLabels(eset)
		if( length(w <- which(!var %in% allv)) ){
			stop("Could not build coef matrix for data entry '", object$key, "':"
					, " phenotypic variable(s) ", str_out(var[w], Inf), " not in the phenotypic data."
					, "\n  [varLabels: ", str_out(allv, Inf), ']')
		}
		
		if( length(var) == 1L ){
			# convert variable into factor if necessary
			if( !is.factor(v <- eset[[var]]) ) 
				v <- factor(v)
			
			if( is_NA(pidx <- pureIndexes(object)) )
				stop("Could not build coef matrix for data entry '", object$key, "': no pure sample definition for data entry.")
			
			t <- droplevels(v[pidx(eset)])
			n <- levels(t)
			
			# check levels
			if( length(w <- which(!n %in% allv)) ){
				stop("Could not build coef matrix for data entry '", object$key, "':"
						, " some levels of phenotypic variable '", var, "'"
						, " [", str_out(n[w], Inf), "] are not in the phenotypic data."
						, "\n  [varLabels: ", str_out(allv, Inf), ']')
			}
		}else n <- var
		
		t(pData(eset)[, n]) / scale			
	}
}

#' \code{pureCoef} builds a coefficient matrix for a data entry where all samples are pure, 
#' using one of its phenotypic -- factor -- variable whose levels define the cell type of origin of each sample.
#' 
#' @rdname GEDdata-internals
pureCoef <- function(var){
	function(object){
		# load eset
		eset <- eset(object)
		
		# check variable
		allv <- varLabels(eset)
		if( length(w <- which(!var %in% allv)) ){
			stop("Could not build coef matrix for data entry '", object$key, "':"
					, " phenotypic variable ", str_out(var[w], Inf), " not in the phenotypic data."
					, "\n  [varLabels: ", str_out(allv, Inf), ']')
		}
		
		# convert as a factor 
		if( !is.factor(v <- eset[[var]]) ) 
			v <- factor(v)
		
		# build base diagonal proportion matrix
		n <- nlevels(v)
		bas <- diag(1, n)
		p <- NA_Matrix(n, ncol(eset), dimnames=list(levels(v), sampleNames(eset)))
		# duplicate columns within each group
		mapply(function(i, j){
			p[, j] <<- bas[, i]
		}, 1:n, split(1:ncol(p), v))
		
		# return proportion matrix
		p
	}
}

#' \code{gedDataCheck} check if all or a given registered dataset can be loaded.
#' @rdname GEDdata-internals
gedDataCheck <- function(key=NULL, load=TRUE, cache=TRUE, save=FALSE, verbose=3){
	
	if( is.null(key) ) key <- gedData()
	res <- sapply(key, function(x){
				e <- ExpressionMix(x, load=load, cache=cache, save=save, verbose=verbose)
				if( verbose ) cat("\n")
				TRUE
			})
	invisible(res)
}
