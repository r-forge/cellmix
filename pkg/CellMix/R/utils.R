# Utility functions 
# 
# Author: Renaud Gaujoux
# Creation: 15 Nov 2011
###############################################################################

#' @include options.R
NULL

"%||%" <- function(a, b) if (!is.null(a)) a else b

nMissing <- function(x) sum(is.na(x))

col <- function(sep, ...) paste(..., collapse=sep)

rotate <- function(x){
	n <- sqrt(nrow(x))
	apply(x, 2L, function(y) as.vector(matrix(as.vector(y), n, byrow=TRUE)))
}

NA_Matrix <- function(n, m=n, ...){
	matrix(as.numeric(NA), n, m, ...)		
}

# lookup/match a column given a vector of values
matchColumn <- function(x, data, fixed=TRUE, no.match=NA){
	mfun <- if( fixed ) is.element else{
		if( length(x) != 1 ) stop("Invalid pattern: must be a single character string [", length(x), ']')
		grepl
	}
	m <- apply(data, 2L, function(v) any(mfun(x,v)))
	w <- which(m)
	if( length(w) > 0L ) w[1L] else no.match
}

#' Detect Log-transformed Data
#' 
#' \code{is_logscale} tells if some numeric data is in log scale, 
#' e.g., normalized microarray data, using the same heuristic as GEO2R.
#' 
#' The data needs to be of reasonable size and variance for the detection 
#' heuristic to work correctly.
#' 
#' @param x a numeric data object (matrix, vector, ExpressionSet) 
#' 
#' @source \url{www.ncbi.org/geo}
#' 
#' @export
#' @examples
#' 
#' x <- rmatrix(20, 10, dist=rnorm, mean=500)
#' is_logscale(x)
#' is_logscale(log_transform(x))
#' 
is_logscale <- function(x){
	
	ex <- if( is(x, 'ExpressionSet') ) exprs(x) else x
	# check log2 transform
	#ex <- exprs(gset)
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	LogC <- (qx[5] > 100) ||
			(qx[6]-qx[1] > 50 && qx[2] > 0) ||
			(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	!LogC
#	if (LogC) { ex[which(ex <= 0)] <- NaN
#		exprs(gset) <- log2(ex) }
}

#' \code{log_transform} apply a log transformation to the data.
#' Negative values get assigned the value \code{\link{NaN}}.
#' 
#' @param base log base to use.
#' 
#' @export
#' @rdname is_logscale
log_transform <- function(x, base=2){
	
	ex <- if( is(x, 'ExpressionSet') ) exprs(x) else x
	
	# transform
	ex[which(ex <= 0)] <- NaN
	ex <- log(ex, base)
	
	# return same type of object
	if( is(x, 'ExpressionSet') ) exprs(x) <- ex
	else x <- ex
	x
}
	

#' Generating Names
#' 
#' @param object object on which to add names
#' @param prefix prefix to use for each name
#' 
#' @return an object of the same class as \code{object} 
#' @export
#' @examples
#' 
#' l <- list(1,2,3,4)
#' addNames(l)
#' addNames( setNames(l, letters[1:length(l)]) )
#' 
addNames <- function(object, prefix=NULL){
	if( is.null(names(object)) )
		names(object) <- paste(prefix,seq_along(object), sep='')
	object
}

# taken from GSEABase:::.checkRequired
checkRequired <- function (required, provided) 
{
	idx <- which(!(required %in% provided))
	if (length(idx) > 0) 
		stop("missing required argument(s): '", paste(required[idx], 
						collapse = "', '"), "'")
}


#' Merge two lists
#' 
#' Adds elements from a list to another list, optionally overwriting them if 
#' already present.
#' 
#' @param x list to which elements from \code{y} are added. 
#' @param y list whose elements are added to \code{x}
#' @param overwrite a logical indicating whether the elements of \code{y} 
#' already present in \code{x} should overwrite the latter or not (default).  
#' 
mergeList <- function(x, y, overwrite=FALSE){
	
	# get the names to add
	extra <- 
	if( !overwrite ) which(! names(y) %in% names(x))
	else seq_along(y)

	# add/overwrite x
	sapply(extra, function(i){
		x[[names(y)[i]]] <<- y[[i]]
	})

	# return updated x
	x
}

#' Returns the index of the first match in a reference table
#' 
#' This function provides a modified behaviour of \code{\link{match}} in the 
#' case of multiple matches, where it returns the index of the first match. 
#' 
#' @inheritParams base::match  
#' 
#' @seealso \code{\link{match}}, \code{\link{charmatch}}
match_first <- function(x, table, nomatch = NA_integer_, incomparables = NULL){
	ndup <- !duplicated(table)
	idx <- seq_along(table)[ndup]
	mi <- match(x, table[ndup], nomatch=0L, incomparables = incomparables)
	res <- rep(nomatch, length(x))
	ma.idx <- mi != 0L
	res[ ma.idx ] <- idx[ mi[ma.idx] ]
	res
}

test.match_first <- function(){	
	x <- c('', 'a', 'ab')
	RUnit::checkIdentical(match_first(x, c('', 'a')), as.integer(c(1, 2, NA)))
	RUnit::checkIdentical(match_first(x, c('a', '')), as.integer(c(2, 1, NA)))
	RUnit::checkIdentical(match_first(x, c('c', '', 'a')), as.integer(c(2, 3, NA)))
	RUnit::checkIdentical(match_first(x, c('ab', 'a')), as.integer(c(NA, 2, 1)))
}


#' Partial Sanitizer for xtable
#' 
#' Text sanitizing function to use in argument \code{sanitize.text.function} 
#' of \code{\link{print.xtable}}.
#' This function is used in \pkg{CellMix} to format the content internal registries 
#' for inclusion in vignettes, as they contain citation commands (\\cite{}).
#' 
#' @param skip character string that contain the character not to sanitize.
#' 
#' @return a function that can be used in argument \code{sanitize.text.function} of
#' \code{print.xtable}.
#' 
#' @seealso \code{\link{print.xtable}}  
#' 
#' @export
#' @keywords internal
#' @examples
#' 
#' library(xtable)
#' print(xtable(gedAlgorithmInfo(show=FALSE), citet=TRUE), sanitize.text.function=xsanitize("\\{}$"))
#'  
xsanitize <- function(skip){
	
	.skip <- skip
	function(str) {
		
		ctypes <- c("\\\\", '$', '>', '<', '|', "{", "}", "%", "&", "_", "#", "^", "~")
		csub <- c('SANITIZE.BACKSLASH', "\\$", '$>$', '$<$', '$|$', "\\{", "\\}", "\\%", "\\&", "\\_", "\\#", "\\verb|^|", "\\~{}")
		stopifnot( length(ctypes) == length(csub) )
		# skip some of the caracters
		iskip <- integer()
		if( !is.null(.skip) && nchar(.skip) > 0L ){
			# special treatment for backslash
			if( grepl("\\", .skip, fixed=TRUE) ){
				iskip <- 1L
				.skip <- gsub("\\", "", .skip, fixed=TRUE)
			}
			.skip <- strsplit(.skip,'')[[1]]
			iskip <- c(iskip, match(.skip, ctypes, nomatch=0L))
			iskip <- iskip[iskip!=0L]
		}
	
		# substitute special characters
		result <- str
		isub <- seq_along(ctypes)
		if( length(iskip) > 0L ) isub <- isub[-iskip]
	#	print(ctypes[isub])
	#	print(x)
		for(i in  isub ){
			result <- gsub(ctypes[i], csub[i], result, fixed = i)
		}
		
		## from xtable:::print.xtable
#		result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
#		result <- gsub("$", "\\$", result, fixed = TRUE)
#		result <- gsub(">", "$>$", result, fixed = TRUE)
#		result <- gsub("<", "$<$", result, fixed = TRUE)
#		result <- gsub("|", "$|$", result, fixed = TRUE)
#		result <- gsub("{", "\\{", result, fixed = TRUE)
#		result <- gsub("}", "\\}", result, fixed = TRUE)
#		result <- gsub("%", "\\%", result, fixed = TRUE)
#		result <- gsub("&", "\\&", result, fixed = TRUE)
#		result <- gsub("_", "\\_", result, fixed = TRUE)
#		result <- gsub("#", "\\#", result, fixed = TRUE)
#		result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
#		result <- gsub("~", "\\~{}", result, fixed = TRUE)
		##
		# substitute basckslashes now
		result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", result, fixed = TRUE)
		return(result)
	}
}

# inspired from GSEABase
slotAccessors <- function(class, slots=slotNames(class), skip=NULL, where=topns(), ...) {
	klass <- class
	slots <- setNames(slots, slots)
	slots <- slots[ !slots %in% skip ]
	for (i in seq(along=slots)) {
		vtype <- getSlots(klass)[[ slots[[i]] ]]
		if (vtype == "ScalarCharacter") {
			vtype <- "character"
			value <- quote(mkScalar(value))
		} else {
			value <- quote(value)
		}
		eval(substitute({
							
			# getter
			if (!isGeneric(GENERIC))
				setGeneric(GENERIC,
						function(object) standardGeneric(GENERIC),
						where=WHERE)
			setMethod(GENERIC,
					signature=signature(object=CLASS),
					function(object) slot(object, SLOT),
					where=WHERE)
			
			# setter
			if (!isGeneric(SETTER) )
				setGeneric(SETTER, function(object, ..., value)
							standardGeneric(SETTER),
						where = WHERE)
			if( '...' %in% names(formals(SETTER)) ){
				setReplaceMethod(GENERIC,
						signature=signature(
								object=CLASS,
								value=VTYPE),
						function(object, check=TRUE, value) {
							slot(object, SLOT) <- VALUE
							if( check ) validObject(object)
							object
						},
						where = WHERE)
			}else{
				
				setReplaceMethod(GENERIC,
						signature=signature(
								object=CLASS,
								value=VTYPE),
						function(object, value) {
							slot(object, SLOT) <- VALUE
							object
						},
						where = WHERE)
				
			}
		}, list(CLASS=klass,
				GENERIC=names(slots)[[i]],
				SETTER=paste(names(slots)[[i]], "<-", sep=""),
				SLOT=slots[[i]],
				VTYPE=vtype,
				VALUE=value,
				WHERE=where)))
	}
}

# set multiple slots at a time from a list of values
setSlots <- function(object, values, skip=NULL, all=TRUE){
	
	if( !all )
		values <- values[names(values) %in% slotNames(object)]
	if( !is.null(skip) )
		values <- values[ ! names(values) %in% skip ]
	
	sl <- getSlots(class(object))[names(values)]
	for( s in names(sl) ){
		slot(object, s) <- as(values[[s]], sl[[s]])
	}
	object
}

#' User Queries
#' 
#' This function is an improved version of \code{\link[Biobase]{userQuery}} and ask 
#' the user about some task that needs her intervention to proceed, 
#' e.g., ask if one should perform a computation, install a package, etc.. 
#' 
#' @inheritParams Biobase::userQuery
#' @param idefault default response in interactive mode.
#' This answer will be in upper case in the question and will be the one returned if the 
#' user simply hits return.
#' 
#' @return the character string typed/agreed by the user or directly the default answer in 
#' non-interactive mode.  
#' 
#' @keywords internal
askUser <- function (msg, allowed = c("y", "n"), idefault = "n", default = "n", case.sensitive = FALSE) 
{
	if ( !interactive())  return(default )
	fallowed <- allowed
	# add extra info on answer options
	if( !is.null(nm <- names(allowed)) ){
		allowed[nm != ''] <- nm[nm != '']  
	}
	if( !isFALSE(idefault) )
		fallowed[fallowed == idefault] <- toupper(idefault)
	repeat {
		allowMsg <- paste("[", paste(fallowed, collapse = "/"), 
				"]: ", sep = "")
		outMsg <- paste(msg, allowMsg)
		cat("\n", outMsg, sep='')
		if (case.sensitive) 
			ans <- readLines(n = 1)
		else ans <- tolower(readLines(n = 1))
		if( !isFALSE(idefault) && !nchar(ans) ) 
			ans <- idefault
		if (ans %in% allowed) 
			break
		else cat(paste(ans, "is not a valid response, try again.\n"))
	}
	# return answer
	ans
}


#' Require a Package with User Interaction
#'
#' Try loading/finding a package and ask the user if it should be installed 
#' if not found.
#' 
#' @param package name of the package
#' @param lib path to the directory (library) where the package is to be
#' looked for and installed if agreed by the user.
#' @param ... extra arguments passed to \code{\link{install.packages}}.
#' @param load a logical that indicates if the package should be loaded, 
#' possibly after installation.
#' @param msg message to display in case the package is not found when first 
#' trying to load/find it.
#' This message is appended to the string \dQuote{Package '<packagename>' is required}.
#' @param quiet logical that indicates if loading a package should be done quietly 
#' with \code{\link{require.quiet}} or normally with \code{\link{require}}.
#' @param prependLF logical that indicates if the message should start at a new line.
#' @param ptype type of package: from CRAN-like repositories, Bioconductor, Bioconductor software, Bioconductor annotation.
#' Bioconductor packages are installed using \code{\link{biocLite}}
#' 
#' @return \code{TRUE} if the package was successfully loaded/found (installed), 
#' \code{FALSE} otherwise.  
#'  
#' @import BiocInstaller
#' @keywords internal
uq_requirePackage <- function(package, lib=NULL, ..., load=TRUE, msg=NULL, quiet=TRUE, prependLF=FALSE
							, ptype=c('CRAN-like', 'BioC', 'BioCSoft', 'BioCann')){
	
	.reqpkg <- if( quiet ) require.quiet else{
				if( prependLF ) message()
				require
			}
	reqpkg <- function(...){
		.reqpkg(..., lib=lib, character.only=TRUE)
	}
	
	# try loading it
	if( load && reqpkg(package) ) return( TRUE )
	# try finding it without loading
	else if( length(find.package(package, lib.loc=lib, quiet=TRUE)) ) return( TRUE )
	
	# package was not found: ask to install
	msg <- str_c("Package '", package, "' is required",
			if( is.null(msg) ) '.' else msg)
	
	# stop if not interactive
	if( !interactive() ) stop(msg)
	
	# detect annotation packages
	if( missing(ptype) && is.annpkg(package) ) ptype <- 'BioCann'
	ptype <- match.arg(ptype)
	
	msg <- str_c(msg, "\nDo you want to install it from known repositories [", ptype, "]?\n"
					, " Package(s) will be installed in '", if(is.null(lib) ) .libPaths()[1L] else lib, "'")
	if( quiet && prependLF ) message()
	repeat{
		ans <- askUser(msg, allowed = c('y', 'n', r='(r)etry'), idefault='y')
		if( ans == 'n' ) return( FALSE )
		if( ans == 'y' ) break
		if( ans == 'r' && reqpkg(package) ) return(TRUE)
	}
	
	## install
	# check Bioconductor repositories
	hasRepo <- function(p){ any(grepl(p, getOption('repos'))) } 
	if( ptype == 'CRAN-like' 
		|| ( ptype == 'BioC' && hasRepo('/((bioc)|(data/annotation))/?$') )
		|| ( ptype == 'BioCsoft' && hasRepo('/bioc/?$') )
		|| ( ptype == 'BioCann' && hasRepo('/data/annotation/?$') )
		){
		install_type <- 'CRAN'
	}
	
	if( install_type == 'CRAN' ){
		pkginstall <- install.packages
	}else{ # Bioconductor 
		if( !reqpkg('BiocInstaller') ){ # get biocLite from bioconductor.org
			# use internal sourceURL to avoid issues with proxies
			sourceURL("http://www.bioconductor.org/biocLite.R")
		}
		# use biocLite wrapper to disable (auto-)updates
		pkginstall <- function(pkgs, ...){
			biocLite(pkgs, ..., suppressUpdates=TRUE, suppressAutoUpdate=TRUE)
		}
	}
	pkginstall(package, lib=lib, ...)
	#
	
	# try  reloading
	if( load ) reqpkg(package)
	else length(find.package(package, lib.loc=lib, quiet=TRUE))
}

# internal source function to play well with CNTLM proxies
sourceURL <- function(url){
	
	file <- url
	if( grepl("^http", url) ){
		dest <- tempfile(basename(url), fileext='.R')
		download.file(url, dest)
		if( file.exists(dest) ){
			file <- dest
			on.exit( file.remove(file) )
		}else stop("Failed to download file '", url, "'")
	}
	source(file)
}

# caching feature adapted from memoise 0.1
new_cache <- function (){
	cache <- NULL
	cache_reset <- function(key) {
		
		if( missing(key) ) cache <<- new.env(TRUE, parent=emptyenv())
		else if( cache$cache_has_key(key) ){
			remove(list = key, envir = cache)
		}
	}
	cache_set <- function(key, value) {
		assign(key, value, envir = cache)
	}
	cache_get <- function(key) {
		get(key, envir = cache, inherits = FALSE)
	}
	cache_has_key <- function(key) {
		exists(key, envir = cache, inherits = FALSE)
	}
	cache_compute <- function(key, value){
		if( cache_has_key(key) ) cache_get(key)
		else{
			value <- force(value)
			cache_set(key, value)
			value
		}
	}
	cache_reset()
	list(reset = cache_reset, set = cache_set, get = cache_get 
		, has_key = cache_has_key, keys = function() ls(cache)
		, compute = cache_compute
		)
}

.cacheXspecies <- new_cache()

add_suffix <- function(x, suffix){
	reg <- gsub(".", "\\.", suffix, fixed=TRUE)
	if( !grepl(str_c(reg, "$"), x) ){
		x <- str_c(x, suffix)
	}
	x
}

comparePackageVersion <- function(pkg, version, ...){
	utils::compareVersion(as.character(packageVersion(pkg, ...)), version)
}

silenceF <- function(f, verbose=TRUE){
	
	if( verbose ) f
	else{
		function(...){
			capture.output(suppressMessages(res <- f(...))); 
			res
		}
	}
}

# force using some -- hidden -- functions from pkgmaker
packageReference <- local({
	.cache <- new_cache()
	function(x, ...){
		.cache$compute(x, pkgmaker:::packageReference(x, ...))
	}
})

# get S3 method even when class is a vector of length > 1
selectS3method <- function(f, class, optional=FALSE){
	m <- lapply(class, getS3method, f=f, optional=TRUE)
	i <- which(!sapply(m, is.null))
	if( !length(i) ){
		if( !optional )
			stop("Could not find a suitable S3 ", f, " method for any of the classes ", str_out(class, Inf))
		return(NULL)
	}
	m[[i[1L]]]
}

# Simple text iteration counter
iterCount <- function(n=100, i0=0, title='Iterations'){
	# setup
	size_i <- nchar(as.character(n))
	.msg_fmt <- paste0(title, ': %', size_i, 'i', "/", n)
	cat(.msg_size <- sprintf(.msg_fmt, i0))
	.msg_size <- nchar(.msg_size)
	.i <- i0
	
	# increment function
	function(i){
		if( is.null(i) ){
			if( .i < n ) cat("\n")
			return()
		}
		if( missing(i) ) i <- .i
		cat(paste0(rep("\r", .msg_size), sprintf(.msg_fmt, i)))
		if( i == n ) cat("\n")
		.i <<- i+1
	}
}