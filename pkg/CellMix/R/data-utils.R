# Functions to handle Datasets
# 
# Author: Renaud Gaujoux
# Created: Apr 1, 2013
###############################################################################

#' @include GEOquery.R
NULL

#' Downloading Gene Expression Deconvolution Datasets
#'  
#' \code{GEDdownload} downloads gene expression datasets that 
#' for datasets that have an entry in the \pkg{CellMix} package's
#' internal registry. 
#' 
#' @param x access key of a registered dataset
#' @param ... extra arguments passed to \code{\link{download.file}} and/or
#' \code{\link{getGSE}}.
#' @param destdir destination directory for all the downloaded files
#' @param datasource Data source from where to fetch the dataset, e.g., 
#' an online dataset repository lke \emph{GEO} or \emph{ArrayExpress}.
#' If missing, it is determined automatically from the access key, but is honoured if provided.
#' Currently only values \code{'GEO'} and \code{'ArrayExpress'} (or its alias \code{'AE'}) are supported.
#' @param annotation annotation package (a string) to attach to the
#' created ExpressionSet object. 
#' @param verbose logical or numeric that indicate the verbosity level
#' 
#' @return the value returned by \code{\link{getGSE}}. 
#' 
#' @export
#' 
GEDdownload <- function(x, ..., destdir=GEDtmp(), datasource=NULL, annotation=NULL, verbose=TRUE){
	
	if( !is.verbose() ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	
	# get registry entry
	dobj <- gedData(x, error=FALSE)
	
	if( is.GEDdata(dobj) ){
		x <- dobj
		if( grepl("^(https?)|(ftp)", x$url) ){
			file <- file.path(destdir, basename(x$url))
			if( !file.exists(file) ){
				download.file(x$url, destfile=file, ...)
			}else message("# Using cached version: '", file, "'")
			message("# Loading file '", file, "'")
			return( getGSE(filename=file, destdir=destdir, ...) )
		}
		key <- x$key
	}else if( isString(x) && is.null(dobj) ){# not found in the local registry: try online databases
		key <- x
	}else{
		stop("Could not download dataset from the given specification [", class(x), ']')
	}
	
	# check type
	datasource <- DataSource(datasource=datasource)
	if( is.null(datasource) ){
		datasource <- DataSource(key)
		if( grepl("^GSE", key) ) 'GEO' # GEO accession number			
		else if( grepl("^E-((MTAB)|(TABM)|(GEOD)|(MEXP))-", key) ) 'ArrayExpress' # ArrayExpress accession number
		else{
			stop("Could not identify data provider of accession key '", key, "'.")
		}
	}	
		
	vmessage("# Fetching dataset '", key, "' from ", datasource)
	vmessage("# Local cache directory: '", destdir, "'")
	if( datasource == 'GEO' ){# GEO accession number (kind of buggy)
		eset <- getGSE(key, destdir=destdir, ...)
		
	}else if( datasource =='ArrayExpress' ){ # ArrayExpress (even more buggy than GEOquery)
		library(ArrayExpress)
		# try local 
		if( lverbose() <= 1 ){ # show warnings if high verbose level
			opt <- options(warn=-1L)
			on.exit( options(opt) )
		}
		vmessage("# Checking local files ...")
		ae <- try(getAE(key, path=destdir, type='processed', local=TRUE), silent=TRUE)
		if( is(ae, 'try-error') ){
			vmessage("# FAILED")
			vmessage("# Downloading data files ...")
			ae <- getAE(key, path=destdir, type='processed', local=FALSE)
			vmessage("# DONE")
		}else vmessage("# OK")
		# FIXME: ArrayExpress - procset crashes if multiple processed files are present, 
		# even if only one of them is actually needed.
		# Here we remove unnecessary files from the list
		# read sample annotation file
		if( (nfiles <- length(ae$processedFiles)) > 1 ){
			vmessage("# Limiting to necessary expression data files ... ", appendLF=FALSE)
			if( file.exists(sannot <- file.path(destdir, paste0(key, '.sdrf.txt'))) ){
				la <- readLines(sannot)
				fused <- sapply(ae$processedFiles, function(f) sum(grepl(f, la[-1], fixed=TRUE)) > 1)
				if( length(not_used <- which(!fused)) ){
					ae$processedFiles <- ae$processedFiles[-not_used]
				}
				vmessage("OK [", length(ae$processedFiles), "/", nfiles,"]")
			}else vmessage("NO [no sample annotation]")
		}
		#
		vmessage("# Loading/Processing dataset ...")
		cnames = getcolproc(ae)
#		print(ae)
#		print(cnames)
		eset <- procset(ae, cnames[2])
		vmessage("# DONE")
		
		# FIXME: ArrayExpress - procset somteimes adds an extra column full of NAs
		if( all(is.na(exprs(eset)[, ncol(eset)])) ){
			vmessage("# Dropping last invalid column ...")
			eset <- eset[, -ncol(eset)]
			vmessage("OK [", ncol(eset), ']')
		}
		# FIXME: ArrayExpress makes a bad job in loading the sample annotation data
		if( !ncol(pData(eset)) ){
			vmessage("# Fixing sample phenotypic annotation ...")
			pd <- ArrayExpress:::readPhenoData(ae$sdrf, destdir)
			pd <- pData(pd)
			sid <- sampleNames(eset)
			if( is_NA(idc <- matchColumn(sid, pd)) ){
				if( is_Illumina(featureNames(eset)) ){ # try removing X prefix
					sid <- substring(sid, 2)
					idc <- matchColumn(sid, pd)
				}
			}
			if( !is_NA(idc) ){
				i <- match(sid, pd[[idc]])
				if( !any(is.na(i)) ){
					pd <- pd[i, ]
					rownames(pd) <- sid
					pData(eset) <- pd
					vmessage("OK [attached ", length(varLabels(eset)), ' variables]')	
				}else{
					vmessage("FAILED [missing some sample annotations]\n  NOTE: Check that the annotation file '", ae$sdrf, "' contains data for all samples.")
				}
			}else vmessage("FAILED [could not match samples]\n  NOTE: Check that the annotation file '", ae$sdrf, "' contains data for all samples.")
		}
		#
	}
	
	# set annotation if provided
	if( !is.null(annotation) && is(eset, 'ExpressionSet') ){ 
		annotation(eset) <- annotation
	}
	
	# return ExpressionSet object
	eset
} 


### DATA SOURCES
DataSource <- function(x, ...){
	UseMethod('DataSource')	
}

#' @S3method DataSource default
DataSource.default <- function(x, datasource=c('GEO', 'ArrayExpress', 'AE'), quiet=FALSE){
	
	# aliases
	aliases <- c(AE='ArrayExpress')
	
	# check validity and aliases for datasource
	if( !missing(datasource) && !is.null(datasource) ){
		datasource <- match.arg(datasource)
		
		# resolve aliases
		if( datasource %in% names(aliases) ){
			datasource <- setNames(aliases[datasource], NULL)
		}
	}else if( missing(datasource) ){
		if( !missing(x) ) datasource <- NULL # will be infered by key
		else{ # list all data sources
			# remove aliases
			datasource <- datasource[! datasource %in% names(aliases)]
		}
	}
	
	# early exit if no key is provided
	if( missing(x) ) return(datasource)
	
	# infer from key if datasource not already defined at this stage 
	if( is.null(datasource) ){
		datasource <- 
			if( grepl("^GSE", x) ) 'GEO' # GEO accession number			
			else if( grepl("^E-((MTAB)|(TABM)|(GEOD)|(MEXP))-", x) ) 'ArrayExpress' # ArrayExpress accession number
			else if( !quiet ){
				stop("Could not identify data provider of accession key '", x, "'.")
			}else return()
	}
	
	# TODO: define S3 structure for data sources (name, regexp key pattern, url prefix, download method, etc...)
	# return data source
	datasource
}

# Returns the url associated with a given GED data.
GEDurl <- function(x){
	
	# get registry entry
	x <- gedData(x)
	
	if( grepl("^((https?)|(ftp))://", x$url) )
		return(x$url)
	key <- x$key
	# determine data source from accession key
	ds <- DataSource(key)
	
	if( ds == 'GEO' ){# GEO accession number
		paste0('http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=', key) 
	}else if( ds == 'ArrayExpress' ){
		paste0("http://www.ebi.ac.uk/arrayexpress/experiments/", key)
	}else
		stop("Could not identify data provider for accession key '", key, "'.")
}

#' Internal Functions for CellMix Data Registry 
#' 
#' These functions are used internally to work with the gene expression dataset
#' registry within the CellMix package.
#' 
#' \code{GEDsave} saves the primary dataset loaded in the registry on disk. 
#' 
#' @param x GEDdata_entry object or list of GEDdata_entry objects.
#' If missing then all entries are used.  
#' @param force logical to force saving dataset(s)
#' 
#' @rdname GEDsave
#' @keywords internal
GEDsave <- function(x, force = FALSE){
	
	.save <- function(d_entry){
		#print(d_entry)
		p <- esetPath(d_entry)
		if( (force || !file.exists(p)) && !is.null(es <- eset(d_entry, load = FALSE)) ){
			message("Saving expression data from ", d_entry$key ," in cache file '", p, "'")
			saveRDS(es, p)
		}
	}
	
	# save all if no data is passed
	if( missing(x) ) x <- gedDataInfo(FALSE)$get_entries()
	
	if( is.GEDdata(x) ) .save(x)
	else if( is.list(x) ) lapply(x, GEDsave, force=force)
	else stop("invalid argument `x`: must be a list or a GEDdata_entry object.")
	
	invisible()
}

setOldClass('GEDdata_entry')

#' User Data Directory
#'  
#' \code{userData} returns the path to a local directory where package-related user data can be stored.
#' Note that a directory is \strong{always} created if necessary (see details).
#' 
#' If in interactive mode, the user is asked if the directory can be created in his home directory,
#' otherwise, or if the user does not allow the creation in his home, the directory is created 
#' in the current R session's temporary directory.  
#' 
#' @param ... path parts passed to \code{\link{file.path}} to be appended to 
#' the main path.
#' @param create logical that indicates if the directory should be created if it does not exists.
#' @param package name of the package associated with the user data path.
#' It is used to prefix the path, within the user R data directory. 
#' 
#' @seealso \code{\link{tempdir}}
#' @export
userData <- local({
			function(..., create=TRUE, package='base'){
				
				p <- file.path(Sys.getenv('HOME'), 'R-data', package)
				
				# ask the user about creating the directory
				if( create && !file.exists(p) ){
					ans <- askUser(str_c("The ", package, " user data directory '", p, "' doen't exist. Do you want to create it?")
							, idefault='y', default='y')
					if( ans == 'n' ){
						p <- file.path(tempdir(), 'R-data', package)
					}
					if( !file.exists(p) ){
						message("Creating user data directory '", p, "'")
						dir.create(p, recursive=TRUE)
					}
				}
				file.path(p, ...)
			}
		})

#' \code{GEDpath} returns/sets the path to the local directory where CellMix data (e.g., cache) are stored.
#' 
#' @inheritParams userData
#' @param reset new value for the local directory
#' 
#' @rdname GEDsave 
#' @export
GEDpath <- local({
			.registryPath <- NULL
			function(..., create=TRUE, reset=NULL){
				# initialise path at first call
				if( is.null(.registryPath) || !is.null(reset) ){
					p <- if( is.character(reset) ) reset[1L] 
							else userData(create=create, package='CellMix')
					.registryPath <<- p
				}
				file.path(.registryPath, ...)
			}
		})

#' \code{GEDtmp} returns the path to the local directory where downloaded data from GEO 
#' are stored.
#' 
#' @rdname GEDsave
GEDtmp <- function(..., create=TRUE){
	p <- GEDpath('GEO', ...)
	if( !file.exists(p) && create ){
		dir.create(p)
	}
	p
}
