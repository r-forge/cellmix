# Gene ID inference and conversion
# 
# Author: Renaud Gaujoux
###############################################################################

#' @include utils.R
#' @include MarkerList-class.R
#' @import annotate
#' @import AnnotationDbi
NULL

library(pkgmaker)

#' Auxiliary function to wrap around mget annotation
#' @keywords internal
.lookUp <- function(keys, map){
	if( !is.list(map) ){
		AnnotationDbi::mget(keys, map, ifnotfound=NA)
	}else{
		res <- NAmap(keys) #setNames(as.list(rep(NA, length(keys))), keys)
		# early exit on empty map
		if( !length(map) ) return(res)
		mk <- keys[keys %in% mappedkeys(map)]
		res[mk] <- map[mk]
		res
	}
}

# define a direct mapping
setClass('GeneIdentifierMap'
	, contains = 'GeneIdentifierType'
	, representation = representation(
		map = 'list' # list of maps
	)
	, prototype = prototype(
		type = mkScalar('Map')
	)
)

GeneIdentifierMap <- function(x=list()){
	if( is_GeneIdentifierMap(x) ) return(x)
	if( !is.list(x) )
		stop("Invalid map: must be a list.")
	if( length(x) && !is.list(x[[1L]]) ){
		x <- list(x)
	}
		
	annotation <- ''
	if( length(x) ){
		
		.names <- function(i, j=i){
			tfrom <- idtype(names(x[[i]]))
			tto <- idtype(unlist(x[[j]]))
			str_c(tfrom,'-', tto)
		}
		annotation <- .names(1L, length(x))
		if( is.null(names(x)) ){
			names(x)[1L] <- annotation
			if( length(x) > 1L ){
				names(x)[2:length(x)] <- sapply(2:length(x), .names)
			}
		}
	}
	new('GeneIdentifierMap', map=x, annotation=annotation)
}

#' Annotation Tools
#' 
#' @description
#' The \pkg{CellMix} package contains a few utility functions to facilitate 
#' working with Bioconductor annotations, which extends or enhance functions
#' available in packages such as \pkg{annotate}.
#' 
#' \code{is.annpkg} tells if an object is the name of an annotation package.
#' 
#' @param x an R object, either a character string or an annotation object.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # check annotation pkg name
#' is.annpkg('aaa.db')
#' is.annpkg(c('bbb.db', 'ccc.db'))
#' is.annpkg(c('ddd', 'eee.db'))
is.annpkg <- function(x) is.character(x) && length(x)>0L && all(grepl("\\.db$", x))

#' \code{is.anndb} tells if an object is an annotation db object such as 
#' \code{hgu133a.db}.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # check AnnotationDb object
#' library(hgu133a.db)
#' is.anndb(hgu133a.db)
is.anndb <- function(x) is(x, 'AnnotationDb')

#' \code{biocann_mapname} returns the name of a map in an annotation package.
#' 
#' @param annotation names of an annotation package, with \dQuote{.db} 
#' suffix or not.
#' @param map name of a map, e.g., \dQuote{ENTREZID}.
#' @param all logical that indicates if all possible names should be 
#' returned, and only the simple concatenation of the annotation 
#' package's name without \dQuote{.db} and the map name.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # build annotation bimap object names
#' biocann_mapname('hgu133a.db', 'ENTREZID')
#' # '.db' extension is skipped
#' biocann_mapname('hgu133a', 'ENTREZID')
#' # get all possible map names
#' biocann_mapname('hgu133a.db', 'ENTREZID', all=TRUE)
biocann_mapname <- function(annotation, map, all=FALSE){
	
	base <- biocann_pkgname(annotation, noext=TRUE)
	sep <- ''
	if( all ) sep <- c(sep, '2')
	paste(base, sep, map, sep='')
	
}
#' \code{biocann_pkgname} returns the name of an annotation package, formated from character strings
#' or extracted from annotation objects.
#' 
#' @param noext logical that indicates if returned package names should 
#' contain the extension '.db'.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' # annotation package name
#' biocann_pkgname('aa')
#' # extract the package name from annotation objects
#' biocann_pkgname(hgu133a.db)
#' biocann_pkgname(hgu133aENTREZID)
biocann_pkgname <- function(x, noext=FALSE){
	
	# extract name from annotation db object
	if( is.anndb(x) ) x <- x$packageName
	if( is(x, 'ProbeAnnDbBimap' ) ) x <- strsplit(x@objTarget, ' ')[[1]][2]
	
	if( !is.character(x) )
		stop("Invalid argument `x`: character string expected [", class(x), ']')
	
	base <- sub("\\.db$", "", x)
	if( noext ) base 
	else ifelse(nchar(base), paste(base, ".db", sep=''), '')
	
}
#' \code{biocann_pkgobject} retrieve the \code{AnnotationDb} object for an annotation
#' package given as a string.
#' The package does not need to be already loaded but needs to be installed in a library that 
#' is in the search path.
#' 
#' @rdname annotation-tools
#' @export
#' @examples
#' 
#' # get AnnotationDb object
#' biocann_pkgobject('hgu133a') # extension can be skipped
#' # the package needs to be installed
#' try( biocann_pkgobject('aaa') )
#' 
biocann_pkgobject <- function(x){
	
	pkg <- biocann_pkgname(x)
	get(pkg, envir=asNamespace(pkg))
	
}

setGeneric('revmap', package='AnnotationDbi')
#' The method \code{revmap} for character vectors revert a mapping provided as 
#' a named character vector.
#' 
#' @param simplify a logical that indicates if the reverse mapping should 
#' be returned as a character vector whenever possible (\code{TRUE}), 
#' i.e. when this mapping is 1:1 -- which can happen even if the original 
#' mapping is not 1:1.
#' Otherwise the mapping is returned as a list.
#' @param ... extra arguments passed to \code{\link{revmap,list-method}}.
#' 
#' @rdname annotation-tools
#' @export
setMethod('revmap', 'character'
		, function(x, ..., simplify=TRUE){
			# check validity of mapping
			if( is.null(names(x)) || (pb <- anyMissing(names(x))) || any(pb <- names(x)=='') )
				stop('revmap - Invalid character map: all elements must have a valid name [', str_out(x[pb], use.names=TRUE), ']')
			
			# convert vector into a list
			x <- 
			if( anyDuplicated(names(x)) ) split(x, factor(names(x), levels=unique(names(x))))
			else as.list(x)
			
			res <- revmap(x, ...)
			# check if one can simplify
			if( simplify && all(sapply(res, length)==1L) ) res <- unlist2(res)
			
			# return result
			res
		}
)

#' Retrieving Bioconductor Annotation Maps
#' 
#' The function \code{biocann_object} retrieves annotation objects, like bimaps, from 
#' Bioconductor annotation packages.
#' It is similar to the function \code{getAnnMap} from the \pkg{annotate} package, 
#' except that it also accepts annotation -- bimap -- objects, 
#' and will try to install missing packages if not found (see section Details).
#' 
#' If an annotation package is specified as a character string, and is not found in the 
#' search path, and if R runs in interactive mode, then the user is asked whether one 
#' should try install the missing package.
#' Default response is 'yes' so that simply hitting return will install the package
#' via \code{\link{install.packages}} and load it.
#' An error is thrown if this eventually fails.
#' 
#' @param to target annotation field as a character string, e.g., \dQuote{ENTREZID}, 
#' \dQuote{ENSEMBL}, or an annotation package or db which means that one wants to 
#' retrieve a mapping to its corresponding primary identifier. 
#' If \code{from} is missing, \code{to} must be the name of an annotation package, 
#' i.e. ends with \dQuote{db}), in which case it tries loading the package and return 
#' the whole annotation db object; or any annotation package db or map object such as 
#' \code{AnnDbBimap}, \code{ChipDb} or \code{OrgDb} objects, which are returned unchanged.
#' @param from source annotation package as a character string e.g. \code{"hgu133a.db"}.
#' @param optional logical that indicates if the function should return \code{NULL} if the 
#' mapping cannot be found (\code{TRUE}), or throw an error.
#' Note that this does not apply to the installation part: if a required annotation 
#' package is missing, an error is thrown even if \code{optional=TRUE}. 
#' 
#' @return a \code{\link{ProbeAnnDbBimap}} if annotation is not missing, 
#' a \code{ProbeAnnDb} object otherwise.
#' 
#' @export
#' @examples
#' 
#' # db package object
#' biocann_object('hgu133a.db')
#' 
#' # bimap from hgu133a probe id to ENTREZID 
#' biocann_object('ENTREZID', 'hgu133a.db')
#' 
#' # reversed bimap from UNIGENE to hgu133a probe id
#' biocann_object('hgu133a.db', 'UNIGENE')
#' # this is equivalent to using the annotation package object (no quotes),
#' # when the package is already loaded (=> helpful in interactive session with auto-completion)
#' biocann_object(hgu133a.db, 'UNIGENE')
#' 
biocann_object <- function(to, from=NULL, optional=FALSE){
	
	if( is.null(from) ){
		# simply return if the argument is already an annotation Bimap
		if( is(to, 'AnnDbBimap') 
				|| is(to, 'ChipDb') 
				|| is(to, 'OrgDb') ){
			return(to)
		}
		# to must be the name of an annotation package 
		if( !is.annpkg(to) )
			stop("Invalid annotation package name ['",to,"']")
		
		# try load and ask if not possible 
		if( !uq_requirePackage(to, load=TRUE) )
			stop("Aborted computation due to missing annotation package '", to, "'.")
		
		ns <- asNamespace(to)
		return(get(to, ns))
		
	}
	
	if( is.character(from) ){
		
		if( is.annpkg(to) || is.anndb(to) ){ # e.g., we want 'ENTREZID -> pkg.db' 
			if( is.annpkg(from) )
				stop("Can only map between an annotation package and a data field.",
						" [map='", to, "' - annotation=", class(from), "]")
			
			# get map and reverse it
			m <- biocann_object(from, to, optional = optional)
			if( !is.null(m) ) return( revmap(m) )
			else return()
		}
		
		# load db package
		annenv <- biocann_object(from)
		
	}else if( is.environment(from) ){
		annenv <- from
	}
	
	if( !is.character(to) )
		stop("Invalid argument `to`: expected character string [", class(to), '] when `from` is [', class(from), ']')
	
	# get all potential map names 
	maps <- biocann_mapname(from, to, all=TRUE)
	for( mname in maps ){
		if( exists(mname, annenv) )	return( get(mname, annenv) )
	}
	
	# error if not optional
	if( !optional ){
		annpkg <- biocann_pkgname(from)
		stop("Could not find map for '", to, "' in package '", annpkg, "'")
	}
	NULL
}

# from GSEABase:::.mapIdentifiers_normalize (1.18.0) 
.mapIdentifiers_normalize <- function (from, to, default=NULL)
{
	from <- asGeneIdentifierType(from)
	to <- asGeneIdentifierType(to)
	
	# early exit if direct mapping
	if( is_GeneIdentifierMap(to) ){
		return( list(from, to) )
	}
	
	#.mapIdentifiers_isMappable(from, to)
	if (nchar(annotation(from)) == 0) 
		annotation(from) <- annotation(to)
	else if (nchar(annotation(to)) == 0) 
		annotation(to) <- annotation(from)
	# use default if possible and necessary
	if( nchar(annotation(from)) == 0  && !is.null(default) ){
		annotation(to) <- annotation(from) <- default
	}
	# add '.db' suffix
	if (nchar(annotation(from)) > 0)
		annotation(from) <- add_suffix(annotation(from), '.db')
	if (nchar(annotation(to)) > 0)
		annotation(to) <- add_suffix(annotation(to), '.db') 
	# return normalized result
	list(from, to)
}

# from GSEABase:::.GeneIdType_asString (1.18.0)
.GeneIdType_asString <- function (x, y) 
{
	if( missing(y) ){
		x <- asGeneIdentifierType(x)
		paste(geneIdType(x), if (nchar(annotation(x)) > 0) {
					paste(" (", annotation(x), ")", sep = "")
				}, sep = "")
	}else{
		sy <- .GeneIdType_asString(y)
		if( is_GeneIdentifierMap(y) ) paste('using', sy)
		else paste('from', .GeneIdType_asString(x), 'to', sy)
	}
}


orgkey <- function(x){
	x <- strsplit(x, " ")[[1]]
	toupper(str_c(substr(x[1],1,3), substr(x[2],1,2)))
}
inp.pkg <- function(x){
	x <- strsplit(x, " ")[[1]]
	p <- str_c(toupper(substr(x[1],1,1)), tolower(substr(x[2],1,1)))
	p <- str_c('hom.', p, '.inp.db')
}

isOrg <- function(x) {
	ann <- annotation(x)
	!is.null(ann) && length(grep("^org\\.", ann)) == 1
}

crossspecies_map <- local({
			
			#.cache <- new_cache()
			.cache <- .cacheXspecies
			function(ez, org1, org2, ifnotfound=NULL, cache=TRUE, lib=NULL, verbose=FALSE){
				
				# set verbosity level
				if( !missing(verbose) ){
					ol <- lverbose(verbose)
					on.exit( lverbose(ol) )
				}
				verbose <- lverbose()
				
				# early exit if empty query
				if( !length(ez) ){
					return(list())
				}
				
				# convert across organism using Bioconductor dedicated function
				inp <- inp.pkg(org1)
				msg <- str_c(' in order to convert IDs from ', org1, ' to ', org2)
				if( !uq_requirePackage(inp, load=TRUE, msg=msg, lib=lib, prependLF=verbose) ){
					vmessage('ERROR');
					stop("Required package '", inp, "' not found.")
				}
				vmessage("# Converting from ", org1, " to ", org2, " ... ", appendLF=FALSE)
				hkey <- digest(list(ez, org1, org2))
				if( cache && .cache$has_key(hkey) ){
					ez2 <- .cache$get(hkey)
					vmessage('OK [cached]')
					res <- ez2
				}else{
					err <- try({
								ez2 <- idConverter(ez
										, srcSpecies=orgkey(org1), destSpecies=orgkey(org2),
										, srcIDType='EG', destIDType="EG")
							}, silent=TRUE)
					if( is(err, 'try-error') ){
						vmessage("ERROR")
						warning("An error occured when converting ids cross-species from ", org1, " to ", org2, ":\n"
								, '  ', err)
						return(list())
					}
					vmessage("OK")
					if( cache ) .cache$set(hkey, ez2)
					res <- ez2
				}
				# fill up unfound values if necessary
				if( !is.null(ifnotfound) ){
					res.partial <- res
					res <- as.list(setNames(rep(ifnotfound, length(ez)), ez))
					res[names(res.partial)] <- res.partial
					res
				}
				#
				res
			}
		})
		

# creates an empty map with given keys
NAmap <- function(x){
	as.list(setNames(rep(NA, length(x)), x))
}

#' Retrieving ID Conversion Maps
#' 
#' Computes a complete list of maps to apply in order to convert from one type 
#' of gene identifier to another.
#' 
#' This function is an extension/enhancement of the function 
#' \code{GSEABase:::.mapIdentifiers_selectMaps}, 
#' which cope with the cases of cross-platform and cross-species conversions.
#' 
#' In particular, when called in an interactive session, this function will notifies 
#' the user about missing required annotation packages, and offers to install them 
#' at runtime.
#' 
#' @param from the source tyep of gene identifier specified as a 
#' \code{\linkS4class{GeneIdentifierType}} object, or a character string, e.g. \code{'ENTREZID'}, 
#' if \code{to} already contains annotation data.
#' @param to the source tyep of gene identifier specified as a 
#' \code{\linkS4class{GeneIdentifierType}} object, or a character string, e.g. \code{'UNIGENE'},
#' if \code{from} already contains annotation data.
#' @param keys character vector of gene identifiers to map. This is only used when performing 
#' cross-species mapping.
#' @param cache logical that indicates if -- lengthy -- cross-species mappings should be cached, 
#' so that subsequent calls with the same set of keys returns much faster.
#' @param link the gene identifier type to use to link between platforms.
#' We \strong{strongly} recommend to leave this argument as default \code{'ENTREZID'}.
#' @param lib path to the library (ie. directory) where to first look and install any required 
#' annotation packages. 
#' This library is temporarly prepended to the .libPath, so that user and R standard library 
#' locations are still looked up.
#' @param verbose logical or integer that indicates the level of verbosity. 
#' 
#' @export
#' @examples 
#' 
#' \dontrun{
#' # intra-platform
#' m <- biocann_map('ENTREZID', AnnotationIdentifier('hgu133a.db'))
#' # cross platform
#' m <- biocann_map(UnigeneIdentifier('hgu133a.db'), EntrezIdentifier('hgu133b.db'))
#' # cross platform, cross-species
#' m <- biocann_map(EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('mogene10stprobeset.db'), keys=as.character(1:1000))
#' }
biocann_map <- function(from, to, keys=NULL, cache=TRUE, link='ENTREZID', lib=NULL, verbose=FALSE){
	
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	verbose <- lverbose()
	# logging functions
	log <- getLogger()
	vmessage <- log$message
	lmessage <- log$lmessage
	#
	
	# add library to path if necessary
	if( !is.null(lib) ){
		olib <- .libPaths()
		on.exit( .libPaths(olib), add=TRUE )
		.libPaths(c(lib, olib))
	}
	#
	
	# format message string (from ... to ...)
	str.fromto <- .GeneIdType_asString(x=from, y=to)
	
	# only keys of `from`
	if( missing(to) ){
		# cannot map if one annotation is missing
		if( !hasAnnotation(from) ){
			vmessage('ERROR');
			stop("Cannot load keys for ", str.fromto, ": no annotation specification.")
		}
		afrom <- getAnnotation(from)
		msg <- str_c(' in order to load keys for ', str.fromto)
		if( !uq_requirePackage(afrom, load=TRUE, msg=msg, prependLF=verbose) ){
			vmessage('ERROR');
			stop("Required package '", afrom, "' not found.")
		}
		# return keys only
		return( keys(biocann_pkgobject(afrom)) )
	}

	# normalize to ensure both have annotations
	type <- .mapIdentifiers_normalize(from, to)
	from <- type[[1]]; to <- type[[2]]
	link <- asGeneIdentifierType(link)
	
	# nothing to do if the types are identical
	if( identical(from, to) ) return()
	
	vmessage("# Loading map(s) ", str.fromto, " ", appendLF=FALSE)
	
	# cannot map if one annotation is missing
	if( !hasAnnotation(from) || !hasAnnotation(to) ){
		vmessage('ERROR');
		stop("Cannot map ", str.fromto, " without annotation specification.")
	} 
	
	# extract annotation packages
	afrom <- getAnnotation(from)
	ato <- getAnnotation(to)
	
	# load required annotation packages
	msg <- str_c(' in order to map IDs ', str.fromto)
	if( !uq_requirePackage(afrom, load=TRUE, msg=msg, prependLF=verbose) ){
		vmessage('ERROR');
		stop("Required package '", afrom, "' not found.")
	}
	if( !uq_requirePackage(ato, load=TRUE, msg=msg, prependLF=verbose) ){
		vmessage('ERROR');
		stop("Required package '", ato, "' not found.")
	}
	#
	
	if( identical(afrom, ato) ){ # intra-platform conversion
		if( isAnnotationIdentifier(from) && isAnnotationIdentifier(to) ){
			# if same annotation package exit with a warning
			if( identical(afrom, ato) ){
				vmessage("WARNING");
				warning("Tried to build a map between 2 identical annotation packages ["
						, .GeneIdType_asString(from),"]")
				vmessage('SKIP')
				return()
			}	
		}
		
		# special case for organism packages and their primary key
		if( is(from, "EntrezIdentifier") && isAnnotationIdentifier(to) && isOrg(to) ){
			k <- keys(biocann_pkgobject(ato))
			res <- list()
			res[[str.fromto]] <- as.list(setNames(k, k))
			res
		}else if( is(to, "EntrezIdentifier") && isAnnotationIdentifier(from) && isOrg(from) ){
			k <- keys(biocann_pkgobject(afrom))
			res <- list()
			res[[str.fromto]] <- as.list(setNames(k, k))
			res
		}else{
			# use function from GSEABase
			res <- GSEABase:::.mapIdentifiers_selectMaps(from, to)
			# rename
			if( is.null(res$second) ){
				res <- setNames(res[1], .GeneIdType_asString(from, to))
			}else{
				names(res) <- c(.GeneIdType_asString(from, AnnotationIdentifier(afrom))
						, .GeneIdType_asString(AnnotationIdentifier(ato), to))
			}
			res
		}
		
	}else{ # cross-platform mapping
		vmessage("[x-platform", appendLF=FALSE)
		# extract organism
		orgfrom <- getAnnMap('ORGANISM', afrom)
		orgto <- getAnnMap('ORGANISM', ato)
		
		if( identical(tolower(str_trim(orgfrom)), tolower(str_trim(orgto))) ){ # same organism
			if( isAnnotationIdentifier(from) && isAnnotationIdentifier(to) ){
				vmessage("/x-probe id] ... ", appendLF=FALSE)
				# conversion between probe ids
				res <- list()
				res[[.GeneIdType_asString(from, link)]] <- biocann_map(from, link, verbose = verbose-1)
				res[[.GeneIdType_asString(link, to)]] <- biocann_map(link, to, verbose = verbose-1)
			}else{
				# conversion via the link identifier
				vmessage("/x-id] ... ", appendLF=FALSE)
				res <- list()
				res[[.GeneIdType_asString(from, link)]] <- biocann_map(from, link, verbose = verbose-1) 
				res[[.GeneIdType_asString(link, AnnotationIdentifier(ato))]] <- biocann_map(link, AnnotationIdentifier(ato), verbose = verbose-1)
				res[[.GeneIdType_asString(AnnotationIdentifier(ato), to)]] <- biocann_map(AnnotationIdentifier(ato), to, verbose = verbose-1)
				res
			}
			
		}else{ # cross-species
			vmessage("/x-species] ... ", appendLF=FALSE)
			if( is.null(keys) ){
				vmessage('ERROR')
				stop("Keys are required for converting IDs across species.")
			}
			res <- list()
			# map source ID to Entrez IDs
			if( !is(from, 'EntrezIdentifier' ) ){
				ez <- biocann_map(from, EntrezIdentifier(), verbose = verbose-1)[[1]]
				if( is.null(ez) )
					stop('Unexpected error: returned null map ', .GeneIdType_asString(from, EntrezIdentifier(afrom)))
				#
				mez <- .lookUp(keys, ez)
				mez <- mappedkeysOnly(mez)
				if( length(mez) ){
					# limit map to strict 1-1 injective mappings
					lmessage(2, "  # Limiting to strict 1-1 mapping ... ", appendLF=FALSE)
					N <- length(mez)
					mez <- idFilterInjective(mez)
					if( verbose >= 2 ) mt <- .mtype(mez, N)
#					print(ez)
					ez <- ez[names(mez)]
					lmessage(2, 'OK [', mt, ']')
					mez <- unique(unlist2(mez))
				}else mez <- character()
				# store map
				res[[.GeneIdType_asString(from, EntrezIdentifier(afrom))]] <- ez
			}else{ # already EntrezIdentifiers
				mez <- keys
			}
			
			# early exit if no keys
			# convert to org2's Entrez IDs (using idconverter)
			ez2 <- crossspecies_map(mez, orgfrom, orgto, ifnotfound=NA, cache=cache, verbose = verbose-2)
			res[[.GeneIdType_asString(EntrezIdentifier(afrom), EntrezIdentifier(ato))]] <- ez2
			# convert from Entrez IDs to destination ID
			if( !is(to, 'EntrezIdentifier' ) ){
				ez2to <- biocann_map(EntrezIdentifier(), to, verbose = verbose-1)
				res[[.GeneIdType_asString(EntrezIdentifier(ato), EntrezIdentifier(ato))]] <- ez2to
			}
			#
			res
			
		}
		
	}
	
	# count number of steps
	nsteps <- length(res) 
	if( length(res) > 1L ){
#		fres <- list()
		.flatten <- function(x){
			mapply(function(val, n){
				if( !is_mapIDList(val) ){
					1
					#fres[[n]] <<- val
				}else sum(.flatten(val))
			}, x, names(x))
		}
		nsteps <- sum(.flatten(res))
#		res <- fres
	}
	
	vmessage('OK [', nsteps, ' step(s)]');
	as.mapIDList(res)
}

as.mapIDList <- function(x){
	class(x) <- c('mapID_List', class(x))
	x
}

is_mapIDList <- function(x){
	is(x, 'mapID_List')
}
is_GeneIdentifierMap <- function(x){
	is(x, 'GeneIdentifierMap')
}


is.multiple <- function(x){
	sapply(x, length) > 1L
}

multkeys <- function(x){
	x[is.multiple(x)]
}

mappedkeysOnly <- function(x){
	x[!is.na(x)]
}

########################
# FILTER STRATEGIES
########################
# setup variable environment within namespace
filterStrategy <- simpleRegistry('filterStrategy', verbose=TRUE)

#' Gene Identifier Filtering Strategies
#' 
#' \code{idFilter} provides access to a small registry that contains
#' a set of filtering strategy functions, which can be used in calls to 
#' the \code{\link{mapIDs}} function.
#' 
#' Filter strategies are functions that take a mapping as its first argument 
#' and returns it after performing some filtering.
#' They should always also have argument \code{...}, so that they can be combined 
#' in sequences.
#' Filter strategies may optionally have an argument \code{.last}, whose value will
#' be filled by \code{mapIDs}: \code{TRUE} if and only if the mapping constitutes 
#' the last step in the identifier conversion.
#' 
#' The following filter strategies are currently defined:
#' 
#' @param name name of the strategy
#' @param ... in \code{idFilter}, these are extra arguments that are used to 
#' pre-build a call to the strategy, which is wrapped into a function.
#' These arguments are only used when \code{.wrap=TRUE}.
#'
#' For the different filtering stratgies, e.g., idFilterAuto, these 
#' arguments are either passed to internal calls to other strategies, 
#' or not used, but required so that filters can be composed. 
#' @param .wrap logical that indicates if the function call should be build 
#' and wrapped into a caller function, which takes only two arguments, 
#' \code{map} and \code{.last}, that will receive the current mapping and 
#' a logical that is \code{TRUE} only when \code{map} is the last mapping
#' in the whole sequence of mappings.
#' 
#' The value of \code{.last} will only be passed to the original filter function 
#' if this one has an argument of the same name.  
#' 
#' @return a function
#' 
#' @export
#' @examples
#' 
#' # list of all available filters
#' idFilter()
#' 
#' # show a given filter function
#' idFilter('1:1')
#' 
#' # composed filter
#' f <- idFilter(c('1:1', 'affy'))
#' 
idFilter <- function(name, ..., .wrap=FALSE){
	
	if( !nargs() ) return( filterStrategy$names() )
#	print(name)
	# simplify name if possible
#	print(name)
	if( is.list(name) && length(name)==1L ) name <- name[[1L]]
	
	if( is.function(name) ) fun <- name 
	else if( is.character(name) ){
		# split composed filters
		name <- str_trim(unlist(strsplit(name, '>')))
		# partial match if necessary
		avail <- filterStrategy$names()
		m <- pmatch(name, avail)
		if( length(i <- which(is.na(m))) )
			stop("Could not find filter strategy '", name[i], "': should be one of ", str_out(avail, Inf))
		name <- avail[m]
		if ( length(name) == 1L ) fun <- filterStrategy$get(name)
		else fun <- lapply(name, idFilter, ..., .wrap=.wrap)
	}else if( is.list(name) ){ # apply idFilter to each element of the list
		
		# if there are names: use name as argument list
		if( !is.null(names(name)) && length(setdiff(names(name), 'name')) > 0L ){
			method <- name
			# gather name specification
			names(method)[names(method) == ''] <- 'name'
			meth <- method[names(method) == 'name']
			method[names(method) == 'name'] <- NULL
			method$name <- meth
			return(do.call('idFilter', c(method, list(...), .wrap=.wrap)))
		}
		fun <- sapply(name, 
				function(x, ...){
					if( is.list(x) ){
						# do not use `...`: all arguments are specified in the element itself 
						do.call(idFilter, c(x, .wrap=.wrap))
					}else idFilter(x, ...) # use `...`: e.g. for list of registered names
				}, ..., .wrap=.wrap)
	}else{
		stop("Invalid argument `name`: must be a character string, a function, or a list of filter specifications [", class(name), ']')
	}

	# create composed filter function that sequentially applies multiple filters
	if( is.list(fun) ){
		if( length(fun) == 1L ){
			fun <- fun[[1L]] # not really a sequence
		}else{
			# return raw list if requested
			if( !.wrap ) return(fun)
			
			.funlist <- fun
			name <- paste(sapply(.funlist, attr, 'name'), collapse=' > ')
			fun <- function(map, ..., .last){
				for(f in .funlist){
					if( f$hasLast ) map <- f(map, .last)
					else map <- f(map)
				}
				map
			}
		}
	}
	
	# ensure the strategy is a function
	if( !is.function(fun) ){
		stop("Invalid filter strategy: function expected [", class(fun), ']')
	}
	if( !.wrap ) return(fun)
	
	# pre-build call
	.cl <- list(quote(fun), quote(map), ...)
	# add .last if necessary
	if( hasLast <- '.last' %in% names(formals(fun)) )
		.cl <- c(.cl, .last=quote(.last))
	.cl <- as.call(.cl)
	#
	
	# return customised function, decorated with the filter's name
	sname <- if( is.function(name) ) '<function>' else name
	ExposeAttribute(function(map, .last){ eval(.cl) }, name=sname, hasLast=hasLast, .VALUE=TRUE, .MODE='r')
}

################
# Strategies
################
#' @details \code{idFilterAll} [key \sQuote{all}]: does not filter at all, and simply 
#' return the input map unchanged.   
#' @param map mapping \code{list} object that maps  source identifiers to another type of
#' identifier.
#' @rdname idFilter
#' @export
idFilterAll <- filterStrategy$set('all', function(map, ...){
	map
})
#' @details \code{idFilterFirstN} [key \sQuote{firstN}]: keeps the first \code{n} matches.
#' @param n maximum number of matches to keep. The default is 1, which only keeps the 
#' first match
#' @param exact.first logical that indicates if exact matches should be prioritized, and 
#' appear as first match.
#'
#' @rdname idFilter
#' @export
idFilterFirstN <- filterStrategy$set('firstN', function(map, n=1L, exact.first=TRUE, ...){
	nl <- sum(sapply(map, length)>n)
	if( nl ){
		log_append('(trunking ', nl, ' maps to 1:', n, ') ')
	}
	
	mapply(function(key, x){
		# look for an exact match and put it in front 
		if( exact.first && !is.na(i <- match(key, x)) ){
			x <- c(x[i], x[-i])
		}
		head(x, n)
	}
	, names(map), map, SIMPLIFY=FALSE)
})
#' @details \code{idFilterOneToOne} [key \sQuote{1:1}]: only keeps single matches.
#'
#' @rdname idFilter
#' @export
idFilterOneToOne <- filterStrategy$set('1:1', function(map, ...){
	if( !length(map) ) return(map) 
	l <- sapply(map, length)
	ok <- l==1L
	if( (n <- length(l) - sum(ok)) ){
		log_append('(dropped ', n, ' 1:2+ maps) ')
	}
	map[ok]
})
#' @details \code{idFilterBiOnetoOne} [key \sQuote{1:1}]: only keeps bidirectional one to one 
#' mapping, i.e. that the reverse mapping is also one to one.
#' In other words it return the largest injective mapping.
#' 
#' @param strict logical that indicates if all intermediate mappings are required to be 
#' \emph{injective} one to one, or only the last mapping step.
#' 
#' @rdname idFilter
#' @export
idFilterInjective <- filterStrategy$set('1-1', function(map, strict=TRUE, ..., .last){
	if( !length(map) ) return(map) 
	# get one-to-one mapping
	res <- idFilterOneToOne(map, ...)
	
	# determine if the mapping should be injective or not
	if( !strict && (!missing(.last) && !.last) ) return(res)
	# apply to revmap as well
	revres <- idFilterOneToOne(revmap(res), ...)
	revmap(revres)
})
#' @details \code{idFilterOneToMany} [key \sQuote{1:*}]: only keeps single matches for 
#' intermediate mappings, but keeps all matches (i.e. single and multiple) when mapping to the 
#' last identifier type.
#'
#' @param .last logical that indicates if \code{map} constitutes the last mapping or an 
#' intermediate mapping.
#'
#' @rdname idFilter
#' @export
idFilterOneToMany <- filterStrategy$set('1:*', function(map, ..., .last){
	if( !length(map) ) return(map)
	if( .last ){
		log_append('(keep all finals) ', appendLF=FALSE)
		map
	}else{
		idFilterOneToOne(map, ..., .last=.last)
	}
})
#' @details \code{idFilterAffy} [key \sQuote{affy}]: enables to filter out secondary affy probes, 
#' i.e., that do not match the primary pattern "^[0-9]+_at$".
#' This filter only applies when the destination identifiers are suitable Affymetrix probes.
#'
#' @param secondary specifies what should be done with secondary probes.
#' If \code{secondary=TRUE} then secondary probes are kept in the mapping but primary probes are prioritized,
#' in the sense that they will appear in front of secondary probes in each element of the mapping.
#' If \code{secondary=FALSE} then secondary probes are removed from the mapping, meaning that identifiers that
#' only mapped to secondary probes are then unmapped.
#' If \code{secondary=NULL} (default), then the filter is applied as when \code{secondary=TRUE}, 
#' \strong{except} when source identifiers are also Affymetrix probes and do not contain any secondary probes,
#' in which case the the filters is applied as when \code{secondary=FALSE}.
#'
#' Any other value for \code{secondary} makes the filter to behave as a pass through, 
#' i.e. the map is returned unchanged.  
#'
#' @rdname idFilter
#' @export
idFilterAffy <- filterStrategy$set('affy', function(map, secondary=NULL, ...){
	
	if( !length(map) ) return(map)
		
	# only apply to Affy chips
	if( idtype(map) != '.Affymetrix' ) return(map)
	# check for presence of secondary probe in the source map
	if( is.null(secondary) ){
		secondary <- 
		if( idtype(names(map)) == '.Affymetrix' ) length(grep("^[0-9]+_at$", names(map), invert=TRUE))>0
		else TRUE
	}
	
	# filter out sprobes if requested
	if( isFALSE(secondary) ){
		k <- unlist2(map)
		k2 <- grep("^[0-9]+_at$", k, value=TRUE)
		diff <- length(k) - length(k2)
		if( diff ) log_append('(dropped ', diff,' 2nd-affy probes) ')
		map <- split(k2, names(k2))
		
	}else if( isTRUE(secondary) ){ 
		# preferably select non s_at probes but do not remove them a priori
		map <- sapply(map, function(x){
			# detect primary probes
			i <- grep("^[0-9]+_at$", x)
			if( length(i) && length(i) != length(x) ){
				x <- c(x[i], x[-i])
			}
			x
		}, simplify=FALSE)
		if( lverbose() ){
			s <- grep("^[0-9]+_at$", unlist2(map), invert=TRUE)
			log_append('(kept ', length(s),' 2nd-affy probes) ')
		}
	}
	else log_append('(kept all affy probes)')

	# return filtered map
	map
})
#' @details \code{idFilterAuto} [key \sQuote{auto}]: this filter is a composed filter, 
#' that successively applies filters \sQuote{affy} and \sQuote{firstN}.
#' It is suitable when one wants to map 1 to 1 as many identifiers as possible.
#'  
#' @rdname idFilter
#' @export
idFilterAuto <- filterStrategy$set('auto', function(map, ..., .last){
	map <- idFilterAffy(map, ..., .last=.last)
	map <- idFilterFirstN(map, ..., .last=.last)
	map
})
#' @details \code{idFilterMAuto} [key \sQuote{mauto}]: this filter is a composed filter, 
#' that successively applies filters \sQuote{affy} and \sQuote{1:*}.
#' @rdname idFilter
#' @export
idFilterMAuto <- filterStrategy$set('mauto', function(map, ..., .last){
	map <- idFilterAffy(map, ..., .last=.last)
	map <- idFilterOneToMany(map, ..., .last=.last)
	map
})
#########################

# returns a string that describe the type of mapping, e.g. 1:1 or 1:1-5
.mtype <- function(k, N=NULL){
	
	s <- NULL
	if( length(k) ){
		if( !is.list(k) ) k <- split(k, names(k))
		# compute map type
		nmult <- sapply(k, length)
		s <- if( sum(nmult>1L) ) str_c('1:', min(nmult), '-', max(nmult), ' = ', sum(nmult)) 
		else '1:1'
	}
	if( !is.null(N) ) s <- str_c(length(k), '/', N, if( !is.null(s) ) str_c(' (', s, ')'))
	s
}

#' Mapping Gene Identifiers
#'
#' @description 
#' This function implements a generic mapping workflow that enables mapping gene 
#' identifers between different types, given a mapper function.
#' The default mapper, \code{\link{biocann_map}}, enables mapping identifiers within or across platform, 
#' as well as across species.
#' 
#' @details
#' The mapper function passed in argument \code{mapper} is responsible for providing a sequence of map(s) 
#' that are sequentially applied, starting with the source gene identifiers in \code{keys}.
#' It must at least have the following 4 arguments:
#' 
#' \describe{
#' \item{from}{source gene identifier type, which will be passed as a \code{\link{GeneIdentifierType}}) object.}
#' \item{to}{source gene identifier type, which will be passed as a \code{\link{GeneIdentifierType}}) object.}
#' \item{keys}{the query keys to map, which may be used to help generate maps specific to the query, e.g., for generating cross-species mappings. 
#' They may also be ignored, if already compiled maps exist, e.g., maps in standard annotation packages,
#' as the returned map do not need to be limited to the query keys.}
#' \item{verbose}{a logical or integer that toggles verbose messages at different levels.}
#' }
#' 
#' @inheritParams convertIDs
#' @param keys character vector of identifiers to map.
#' @param method a specification of the filtering strategy. See \code{\link{idFilter}}.
#' @param mapper mapper function, responsible to generate the actual mappings between
#' identifiers. See details.
#' 
#' @export
#' @examples
#' # some ENTREZ ids
#' ez <- as.character(1:1000)
#'
#' # convert into probe ids for Affymetrix Chip hgu133a
#' m <- mapIDs(ez, EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('hgu133b.db'), verbose=2)
#'
#' # keep primary affy probes only
#' m <- mapIDs(ez, EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('hgu133b.db'), method='affy', verbose=2)
#'
#' # same but only keep 1:1 mapping, using a composed filtering strategy
#' m <- mapIDs(ez, EntrezIdentifier('hgu133a.db'), AnnotationIdentifier('hgu133b.db'), method=c('affy', '1:1'), verbose=2)
#' 
mapIDs <- function(keys, from, to, method=c('auto', 'mauto', 'all', 'firstN', '1:1', '1:*', 'affy'), mapper=biocann_map, verbose=FALSE, ...){
	
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	verbose <- lverbose()
	# logging functions
	log <- getLogger()
	vmessage <- log$message
	lmessage <- log$lmessage
	#
	
	maps <- NULL
	# direct mapping
	if( is.list(to) ) to <- GeneIdentifierMap(to)
	
	if( is_GeneIdentifierMap(to) ){ # `to` contains the list of mapping
		if( missing(from) ) from <- GeneIdentifierMap() 
		maps <- to@map
	}else{
	
		# normalize to ensure both have annotations
		type <- .mapIdentifiers_normalize(from, to)
		from <- type[[1]]; to <- type[[2]]
		# extract annotation packages
		afrom <- getAnnotation(from)
		
	}
	
	vmessage("# Converting ", .GeneIdType_asString(from, to), ' ... ', appendLF=verbose>1)
	
	# load filtering strategy
	filter <- NULL
	if( missing(method) ) method <- NA
	if( is.null(method) ) method <- 'auto'
	if( !is_NA(method) ){
		if( !is.list(method) ){
			method <- list(method)
		}
		if( is.null(names(method)) ) names(method) <- rep('', length(method))
		# force wrapping into pre-compiled filter function
		method$.wrap <- TRUE
		# gather name specification
		names(method)[names(method) == ''] <- 'name'
		meth <- method[names(method) == 'name']
		method[names(method) == 'name'] <- NULL
		method$name <- meth 
		# skip filtering if method is 'all'
		filter <- do.call('idFilter', method)
		filterName <- str_out(filter$name)
	}
	#
	
	# init vector-map
	mkeys <- setNames(keys, keys)
	# remove duplicated keys
	if( anyDuplicated(keys) ){
		lmessage(2, "# Removing duplicates from query ... ", appendLF=FALSE)
		mkeys <- mkeys[!duplicated(mkeys)]
		lmessage(2, '[', length(keys), ' -> ', length(mkeys), ' id(s)]')
	}
	# initialize list-map
	mapping <- NAmap(mkeys)
	
	# limit source ids to those from the source annotation package
	lmessage(2, "# Limiting query to ", .GeneIdType_asString(from), ' ... ', appendLF=FALSE)
	if( !is.null(maps) ){
		if( !length(maps) ){
			warning("Did not map any id ", .GeneIdType_asString(from, to),": mapping is empty.")
			return(mapping)
		}
		mkeys <- mkeys[mkeys %in% names(maps[[1L]])]
	}else if( !isAnnotationIdentifier(from) ){
		annmap <- mapper(from=from, to=AnnotationIdentifier(afrom), verbose=verbose-2, ...)[[1]]
		mkeys <- mkeys[mkeys %in% mappedkeys(annmap)]
	}else{
		annmap_keys <- mapper(from=from, to=, verbose=verbose-2, ...)
		mkeys <- mkeys[mkeys %in% annmap_keys]
	}
	lmessage(2, '[', length(mapping), ' -> ', length(mkeys), ' id(s)]')
	#
	# early exit if there are no keys left
	if( !length(mkeys) ){
		vmessage('SKIP', ' [all ', length(keys), ' keys fail filtering]')
		return(mapping)
	}
	# nothing to do if from==to
	if( identical(from, to) ){
		vmessage('SKIP', ' [kept ', length(mkeys), '/', length(keys),' key(s)]')
		mapping[names(mkeys)] <- mkeys
		return(mapping)
	}
	
	# retrieve full map if necessary
	if( is.null(maps) ){
		maps <- mapper(from=from, to=to, keys=mkeys, verbose=verbose-1, ...)
	}
	
	if( !length(maps) ){
		warning("Did not map any id ", .GeneIdType_asString(from, to),": mapping is empty.")
		return(mapping)
	}
	
	.mapLookUp <- function(keys, map, mapname, islast){
		oi <- log$nindent(1L)
		on.exit(log$nindent(oi, add=FALSE))
		lmessage(2, "# Mapping ", mapname, ' [', length(map), ' entries] ... ',  appendLF=FALSE)
		# early skip if empty query
		if( !length(keys) ){
			vmessage('SKIP [0 keys]')
			return(keys)
		}
		k <- .lookUp(keys, map)
		# reset names using primary source ID
		k <- unlist2(setNames(k, names(keys)))
		# remove unmapped
		k <- mappedkeysOnly(k)
		# reshape into a list
		k <- split(k, names(k))
		# remove duplicated matches
		k <- sapply(k, unique, simplify=FALSE)
		if( lverbose() > 1L ){
			# compute map type
			mtype <- .mtype(k)
			vmessage('[', length(k), '/', length(keys), ' mapped (', mtype, ')]')
		}
		# Apply filtering
		if( !is.null(filter) && length(k) ){
			lmessage(2, "# Applying filtering strategy ", filterName, " ... ", appendLF=FALSE)
			nk <- length(k)
			k <- filter(map=k, .last=islast)
			if( lverbose() > 1L ){
				# check how many passed
				mtype <- .mtype(k)
				vmessage('[', length(k), '/', nk, ' passed (', mtype, ')]')
			}
		}
		#

		# return as a vector map
		if( length(k) ) unlist2(k)
		else character()
	}
	
	# walk through the list of maps
	.map_along <- function(maps, islast=TRUE){
		mapply(function(i, mname){
			map <- maps[[i]]
			lastelmt <- i==length(maps)
			if( !is_mapIDList(map) ){
				mkeys <<- .mapLookUp(mkeys, map, mname, islast && lastelmt)
			}else{
				if( !length(map) || length(map) > 2L ){
					warning("Map list has an unexpected length [", length(map), "]")
				}
				.map_along(map, lastelmt)
			}
		}, seq_along(maps), names(maps))
	}
	.map_along(maps)
	#
	
	# build result mapping
	mapping <- mkeys <- split(mkeys, names(mkeys))
	# extend mapping to full query
#	mapping <- sapply(keys, function(k){
#				m <- mkeys[[k]]
#				if( is.null(m) ) NA else m
#			}, simplify=FALSE)
	if( lverbose() ){
		# compute map type
		mtype <- .mtype(mapping)
		vmessage('OK', ' [', length(mkeys), '/', length(keys), ' mapped (', mtype, ')]')
	}
	mapping
	#
	
}

#' Convert nuID to Nucleotide Sequence
#' 
#' The function \code{nuIDdecode} converts a nuID string back to 
#' the nucleotide sequence it encodes.
#' nuIDs are identifiers used as primary keys in some Illumina annotation packages, 
#' and are based on a hash of the probe sequence itself. 
#' 
#' This function is an adapted version of the \code{lumi::id2seq} from the 
#' \pkg{lumi} package, so that it can throws errors instead of warnings.
#' It is used in \code{\link{idtype}} to infer the type of nuID vectors.
#' 
#' @param id nuID character vector
#' @param error a logical indicating whether an error should be thrown in case 
#' of invalid input nuID (default) or only a warning, as in the original function
#' \code{lumi::id2seq}.
#' 
#' @source function \code{lumi::id2seq} in the \pkg{lumi} package
#' 
#' @author Pan Du
#' 
#' Adaptation for CellMix: Renaud Gaujoux
#' @cite Du2007
#' 
#' @export
#' @examples
#' nuIDdecode('XERDqYYc2A')
#' try(nuIDdecode('XERDqYYc2F'))
#' nuIDdecode('XERDqYYc2F', error=FALSE)
#' 
nuIDdecode <- function (id, error = TRUE) 
{
	if (length(id) > 1) {
		return(sapply(id, nuIDdecode, error=error))
	}
	else {
		code <- c(LETTERS, letters, as.character(0:9), "_", ".")
		ind <- 1:length(code)
		names(ind) <- code
		if (length(grep("[^a-zA-Z0-9_.]", id)) > 0){			
			msg <- "Input id is not a nuID!"
			if( isTRUE(error) ) stop(msg)
			else return(NA) 
			#else warning(msg)
		}
		id <- substring(id, 1:nchar(id), 1:nchar(id))
		num <- as.numeric(ind[id]) - 1
		checkCode <- num[1]
		num <- num[-1]
		if (checkCode == 63){
			msg <- "Coding error or not a nuID!\n Check code should not include \".\"!"
			if( isTRUE(error) ) stop(msg)
			else if( is.na(error) ) return(NA) 
			else warning(msg)
		}
		cutLen <- checkCode%%3
		res <- floor(checkCode/3)
		codon <- rbind(floor(num/4^2), floor((num%%4^2)/4), num%%4)
		checkSum <- sum(num)
		if (res != checkSum%%21) {
			msg <- "Coding error or not a nuID!"
			if( isTRUE(error) ) stop(msg)
			else if( is.na(error) ) return(NA)
			else warning(msg)
		}
		nucleotide <- c("A", "C", "G", "T")
		codon <- nucleotide[codon + 1]
		len <- length(codon)
		seq <- paste(codon[1:(len - cutLen)], collapse = "")
		return(seq)
	}
}

###########
# IDTYPE
###########

#' Method for when \code{idtype} is called with its first argument missing, 
#' in which case it returns all or a subset of the known type names as a character 
#' vector, or optionally as a list that contains their definition, i.e. a regular 
#' expression or a matching function.
#' 
#' @param def a logical or a subsetting vector, used when \code{object} is missing, 
#' which indicates that the result should contain the definition of the matching 
#' pattern/function of each type, or which type's deifnition should be included
#' in the result list.
#'   
#' @examples 
#' 
#' # all known types
#' idtype()
#' # with their definitions
#' idtype(def=TRUE)
#' idtype(def='ENTREZID')
#' idtype(def=c('ENTREZID', 'ENSEMBLTRANS'))
#' 
setMethod('idtype', 'missing' 
		, local({
					
					# type matching patterns/functions
					.defs <- list(
							UNIGENE="^[A-Z][a-z]\\.[0-9]+$"
							, ENSEMBL="^ENSG[0-9]+$"
							, ENSEMBLTRANS="^ENST[0-9]+$"
							, ENSEMBLPROT="^ENSP[0-9]+$"
							, ENTREZID="^[0-9]+$"
							, IMAGE = "^IMAGE:[0-9]+$"
							, GOID="^GO:[0-9]+$"
							, PFAM="^PF[0-9]+$"
							, REFSEQ="^[XYN][MPR]_[0-9]+$"
							, ENZYME="^[0-9]+(\\.(([0-9]+)|-)+){3}$"
							, MAP="^(([0-9]{1,2})|([XY]))((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9]{1,2})|([XY]))?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$"
							, GENEBANK=c("^[A-Z][0-9]{5}$", "^[A-Z]{2}[0-9]{6}$")
							, GENEBANK="^[A-Z]{3}[0-9]{5}$"
							, GENEBANK="^[A-Z]{4}[0-9]{8}[0-9]?[0-9]?$"
							, GENEBANK="^[A-Z]{5}[0-9]{7}$"
							, GENENAME = " "
							, .Affymetrix="(^AFFX[-_])|(^[0-9]+_([abfgilrsx]_)?([as]t)|(i))$"
							, .Illumina="^ILMN_[0-9]+$"
							, .Agilent = "^A_[0-9]+_P[0-9]+$"
							, .nuID=function(x) !is.na(nuIDdecode(x, error=NA))
					#, SYMBOL="^[-0-9a-zA-Z':@.]*[a-zA-Z][-0-9a-zA-Z':@.]*$"
					)
					
					# actual function
					function(object, def=FALSE){
						if( isFALSE(def) ) names(.defs)
						else if( isTRUE(def) ) .defs
						else if( is.vector(def) ){
							if( is.character(def) && length(def) == 1L ){
								res <- .defs[names(.defs) %in% def]
								if( length(res) == 1L ) res[[1L]]
								else res
							}else .defs[def]
						}else
							stop("idtype - Invalid argument `def` of class:'",class(def),"': expecting FALSE, TRUE or a subsetting vector.")
					}
				})
)
#' Detects the type of identifiers used in the row names of a matrix. 
setMethod('idtype', 'matrix', function(object, ...) idtype(rownames(object), ...) )
#' Detects the type of identifiers used in the feature names of an \code{ExpressionSet} object. 
setMethod('idtype', 'ExpressionSet', function(object, ...) idtype(featureNames(object), ...) )
#' Detects the type of identifiers used in the rownames of the basis matrix of 
#' an \code{NMF} model.
setMethod('idtype', 'NMF', function(object, ...) idtype(rownames(object), ...) )
#' Detects the type of the primary identifiers of a probe annotation bimap object.
#'  
#' To speedup the identification, only the first 500 probes are used by default, 
#' since the IDs are very likely to have been curated and to be of the same type.
#' This can be changed using argument \code{limit}.
#' 
setMethod('idtype', 'ProbeAnnDbBimap', function(object, limit=500L, ...) idtype(keys(object), limit=limit, ...) )
#' Detects the type of the identifiers of a chip annotation object.
#' 
#' To speedup the identification, only the first 500 probes are used by default, 
#' since the IDs are very likely to have been curated and to be of the same type.
#' This can be changed using argument \code{limit}.
setMethod('idtype', 'ChipDb', function(object, limit=500L, ...) idtype(keys(object), limit=limit, ...) )
#' Returns the type of identifier defined by a \code{\link{GeneIdentifierType}}
#' object.
#' Note that this methods is a bit special in the sense that it will 
#' return the string \dQuote{ANNOTATION} for annotation based identifiers, 
#' but will not tell which platform it is relative to.
#' This is different to what \code{idtype} would do if applied to the primary 
#' identifiers of the corresponding annotation package.
#' 
#' @examples
#' # from GeneIdentifierType objects
#' idtype(NullIdentifier())
#' idtype(AnnotationIdentifier('hgu133a.db'))
#' # but 
#' \dontrun{
#' 	library(hgu133a.db)
#' 	idtype(hgu133a.db)
#' }
#' 
setMethod('idtype', 'GeneIdentifierType', 
		function(object){
			t <- toupper(geneIdType(object))
			if( t == 'NULL' ) t <- ''
			t
		}
)
#' Detects the type of all elements in a list, but provides the option of 
#' detecting the type of each element separately.
#' 
setMethod('idtype', 'list', function(object, ...) idtype(unlist(object), ...) )
setMethod('idtype', 'MarkerList', function(object, each=FALSE, ...){
			if( !length(object) ) return( '' ) 
			if( each ) sapply(dropvalues(object), idtype, each=each>1L, ..., simplify=NA)
			else idtype(marknames(object), each=FALSE, ...)
		})
#' Dummy method -- defined for convenience -- that returns \code{''}
setMethod('idtype', 'NULL', function(object, ...) '' )


is_Affymetrix <- function(x){ idtype(x) == '.Affymetrix'}
is_Illumina <- function(x){ idtype(x) == '.Illumina'}
is_Agilent <- function(x){ idtype(x) == '.Agilent'}

# Returns the manufacturer
str_platform <- function(object){
	t <- idtype(object)
	if( t == '.nuID' ) t <- '.Illumina'
	p <- sub("^\\.", '', t)
	s <- NULL
	if( t != p ) s <- p
	if( hasAnnotation(object) ){
		if( !is.null(s) ) s <- paste(s, '/', str_c(annotation(object), collapse=', '))
		else s <- annotation(object)
	}
	s
}

#' This is the workhorse method that determine the type of ids contained in 
#' a character vector. 
#' 
#' @param each logical indicating whether the type of each element should be 
#' returned (\code{TRUE}) or only the type of the vector as a whole (default).
#' @param limit specification for limiting which elements are used to 
#' detect the type of identifiers.
#' If a single numeric, then only the first \code{limit} elements 
#' are used. Otherwise it must be a subsetting logical or numeric vector.
#' @param no.match character string that specifies the string to use when the 
#' type cannot be determined.    
#' 
#' The IDs can be either:
#' \itemize{
#' \item{probe IDs (e.g. 123456_at or ILMN_123456 for Affymetrix or Illumina 
#' chips respectively), the type starts with a dot \code{'.'}, allowing the 
#' subsequent handling of such IDs as a group.}
#' \item{other biological ID types, the result are character strings such as 
#' those used as attributes in Bioconductor annotation packages 
#' (e.g. \code{"ENTREZID"} or \code{"ENSEMBL"})}
#' \item{Names of annotation packages e.g. \code{"hgu133plus2.db"}.}
#' }
#' 
#' This function is able to identify the following ID types using regular 
#' expression patterns or dedicated function:
#' \itemize{
#' \item ENSEMBL = "^ENSG[0-9]+$"
#' \item ENSEMBLTRANS = "^ENST[0-9]+$"
#' \item ENSEMBLPROT = "^ENSP[0-9]+$"
#' \item ENTREZID = "^[0-9]+$"
#' \item IMAGE = "^IMAGE:[0-9]+$"
#' \item GOID = "^GO:[0-9]+$"
#' \item PFAM = "^PF[0-9]+$"
#' \item REFSEQ = "^N[MP]_[0-9]+$"
#' \item ENZYME = "^[0-9]+(\\.(([0-9]+)|-)+){3}$"
#' \item MAP = "^[0-9XY]+((([pq])|(cen))(([0-9]+(\\.[0-9]+)?)|(ter))?(-([0-9XY]+)?(([pq]?)|(cen))((ter)|([0-9]+(\\.[0-9]+)?))?)?)?$"
#' \item GENEBANK (Nucleotide) = "^[A-Z][0-9]{5}$" | "^[A-Z]{2}[0-9]{6}$" 
#' \item GENEBANK (Protein) = "^[A-Z]{3}[0-9]{5}$"
#' \item GENEBANK (WGS) = "^[A-Z]{4}[0-9]{8}[0-9]?[0-9]?$"
#' \item GENEBANK (MGA) = "^[A-Z]{5}[0-9]{7}$"
#' \item GENENAME = " "
#' \item .Affymetrix = "(^AFFX-)|(^[0-9]+_([abfgilrsx]_)?([as]t)|(i))$"
#' \item .Illumina = "^ILMN_[0-9]+$"
#' \item .Agilent = "^A_[0-9]+_P[0-9]+$"
#' \item .nuID = use the function \code{\link{nuIDdecode}} to try converting the 
#' ids into nucleotide sequences. Identification is positive if no error is 
#' thrown during the conversion.  
#' }
#' 
#' @return a single character string (possibly empty) if \code{each=FALSE} (default) 
#' or a character vector of the same "length" as \code{object} otherwise.
#' @export
#' @examples
#' 
#' idtype("12345_at")
#' idtype(c("12345_at", "23232_at", "555_x_at"))
#' # mixed types
#' ids <- c("12345_at", "23232_at", "Hs.1213")
#' idtype(ids) # not detected
#' idtype(ids, each=TRUE)
#' 
setMethod('idtype', 'vector',
		function(object, each=FALSE, limit=NULL, no.match=''){
			
			no.match <- as.character(no.match)
			
			if( length(object) == 0 ) return(no.match)
			# handle annotation package names
			if( is.character(object) && all(is.annpkg(object)) ){
				res <- 
						if( length(object) == 1L ) idtype(biocann_object(object))
						else sapply(object, function(x){ idtype(biocann_object(x), no.match=NA) })
				return(res)
			}
			
			ids <- object
			# remove missing
			ids <- ids[!is.na(ids)]
			# early exit if empty
			if( length(ids) == 0 ) return(no.match)
			
			# use only a subset of ids 
			if( !is.null(limit) ){
				if( !is.numeric(limit) && !is.logical(limit) )
					stop("idtype - Invalid value for argument `limit`: must be numeric or logical vector.")
				ids <-
						if( length(limit) == 1L ) head(ids, limit)
						else if( is.logical(limit) ) ids[limit]
						else ids[as.integer(limit)]
			}
			
			if( length(ids) == 0 ){
				warning("idtype - Empty input vector: returned NA.")
				return(as.character(NA))
			}
			patterns <- idtype(def=TRUE)		
			
			# special case for integer vectors
			if( is.integer(ids) ){
				if( each ) return(setNames(rep("_INDEX_", length(ids)), ids))
				else return("_INDEX_")
			}
			
			.check <- function(def, keys){
				if( is.character(def) ){
					res <- grepl(def[1L], keys)
					if( length(def)>1L ){
						for( i in 2:length(def) ){
							res <- res | grepl(def[i], keys)
						}
					}
					res
				} else if( is.function(def) ) def(keys)
				else stop("Invalid ID matcher object [", class(def), ']')
			}
			
			if( each ){
				res <- setNames(rep(NA, length(ids)), ids)
				idx <- seq_along(res)
			}
			for( i in seq_along(patterns) ){
				p <- patterns[[i]]
				
				if( !each ){
					if( all(.check(p, head(ids, 3))) && all(.check(p, ids)) )
						return(names(patterns)[i])
				}else{
					t <- .check(p, ids[idx])
					ok <- which(t)
					if( length(ok) > 0 ){
						res[idx[ok]] <- names(patterns)[i]
						idx <- idx[-ok]
						if( length(idx) == 0L ) break
					}
				}
			}
			
			if( each ) res
			else no.match
		}
)


#' Utility function for Biological Identifiers
#' 
#' \code{is.probeid} tells if given IDs are probe IDs.
#' 
#' @param x an R object
#' @param ... extra arguments passed to \code{\link{idtype}}.
#' 
#' @family biocIDs
#' @rdname biocIDs
#' @export
is.probeid <- function(x, ...) is.probetype(idtype(x, ...))
#' \code{is.probetype} tells if given types are probe ID types.
#' 
#' @param type an identifier type as returned by \code{\link{idtype}}.
#' 
#' @rdname biocIDs
#' @export
is.probetype <- function(type) type == 'ANNOTATION' | (is.idtype(type) & grepl("^\\.", type))
#' \code{is.idtype} tells if a given character vector contains valid types.
#' 
#' @rdname biocIDs
#' @export
is.idtype <- function(type){
	if( !is.character(type) ) FALSE
	else type %in% idtype()
}

#' \code{asGeneIdentifierType} is similar to \code{\link{idtype}}, 
#' but returns an object of class \code{\linkS4class{GeneIdentifierType}}, 
#' for use in \code{\linkS4class{GeneSet}} objects.
#' 
#' @inheritParams idtype
#' @param annotation annotation package to associate with the 
#' returned \code{\linkS4class{GeneIdentifierType}} object.
#' @param error specifies what to do if the type of \code{x} cannot be determined 
#' 
#' For \code{asGeneIdentifierType}, all arguments in \code{...} are used in 
#' an internal call to \code{idtype}.
#' 
#' @rdname biocIDs
#' @export 
asGeneIdentifierType <- function(x, limit=NULL, annotation=NULL, error=NullIdentifier(annotation)){
	
	# known idtypes 
	known <- c('UNIGENE', 'ENSEMBL', 'ENTREZID', 'REFSEQ', 'ENZYME', 'GENENAME', 'SYMBOL')
	if( missing(x) ) return(known)
	
	noAnnotation <- missing(annotation)
	
	# do nothing if x is already of the right type
	if( is(x, 'GeneIdentifierType') ){
		if( !is.null(annotation) ){
			x <- new(class(x), annotation=.glue_anndb(annotation))
		}
		return( x )
	}
	# for a list 
	
	type <- NULL
	if( is.annpkg(x) ){
		# check for direct name of type in `x`
		type <- x
	}else if( isString(x) ){
		# check known types (partial matching)
		t <- unique(idtype())
		if( !is.na(m <- pmatch(x, t)) ){
			type <- x <- t[m]
		}else if( !is.na(m <- pmatch(x, known)) ){
			type <- x <- known[m]
		}
	}else {
		# check for an embedded, non trivial, GeneIdentifier object in `x`
		if( hasMethod('geneIdType', class(x)[1L]) && !identical(idtype <- geneIdType(x), NullIdentifier()) ){
			return(idtype)
		}
	}
	
	# check for embedded annotations and use it if not overwriten by argument `annotation`
	if( missing(annotation) && !is.null(ann <- getAnnotation(x, null=TRUE)) ){
		annotation <- ann
	}
	# try directly inferring from the data
	if( is.null(type) ){
		# special case of ProbeAnnDbBimap objects
		if( is(x, 'ProbeAnnDbBimap' ) ){
			type <- idtype(revmap(x))
			if( missing(annotation) ){ 
				annotation <- biocann_pkgname(x)
			}
		}else type <- idtype(x, limit=limit)
	} 
	# use empty annotation string
	if( is.null(annotation) ) annotation <- ''
	annotation <- .glue_anndb(annotation)
	
	# early exit for annotation packages and probe ids
	if( is.annpkg(type) ) return( AnnotationIdentifier(.glue_anndb(type)) )
	if( all(is.probetype(type)) ) return( AnnotationIdentifier(annotation) )
	
	if( !isString(type) ){
		stop("Invalid inferred id type for `", deparse(substitute(x))
			, "` (annotation: '", annotation, "'): should be a single character string [", class(type), ']')
	}
	
	# try other types (partial matching)
	if( type != '' ){
		m <- pmatch(type, known)
		# not found?
		if( is.na(m) ){
			if( isTRUE(error) ){
				if( identical(type, '') ){
					err <- str_c("Could not inferred id type from `", deparse(substitute(x)), "`  (annotation: '", annotation, "')")	
				}else{
					err <- str_c("Invalid inferred id type '", type, "' from `", deparse(substitute(x)), "`  (annotation: '", annotation, "'):"
								, " should be one of ", str_out(known, Inf))
				}
				stop(err)
			}
			return(error)
		}
		type <- known[m]
	}else if( annotation != '' ){
		return( AnnotationIdentifier(annotation) )
	}
	
	if( type == '' ){
		if( isTRUE(error) ){
			if( isString(x) )
				stop("Could not infer type from string '", x, "'  (annotation: '", annotation, "').\n"
				, "  String input should be an annotation package name (*.db) or one of ", str_out(known, Inf))
			stop("Could not infer type from object of class '", class(x), "'.")
		}
		return(error)
	} 
	
	switch(type,
			UNIGENE ={
				if( noAnnotation && is.character(x) && !identical(x, 'UNIGENE') )
					annotation <- str_c('org.', substr(x[[1L]], 1, 2), '.eg.db')
				UnigeneIdentifier(annotation)
			},
			ENSEMBL = ENSEMBLIdentifier(annotation),
			ENTREZID = EntrezIdentifier(annotation),
			REFSEQ = RefseqIdentifier(annotation),
			ENZYME = EnzymeIdentifier(annotation),
			GENENAME = GenenameIdentifier(annotation),
			SYMBOL = SymbolIdentifier(annotation),
			{
				stop("Invalid inferred type '", type, "' (annotation: '", annotation, "'): should be an annotation package name (*.db) or one of ", str_out(known, Inf))
			}
	)
}

mustConvertIDs <- function(x, y, returnType=FALSE){
	itx <- asGeneIdentifierType(x)
	ity <- asGeneIdentifierType(y)
	must <- !isNullIdentifier(itx) && !isNullIdentifier(ity) && !identical(itx, ity)
	if( !returnType ) must
	else if( must ) ity
	else NULL
}

tryConvertIDs <- function(x, y, ...){
	if( !is.null(it <- mustConvertIDs(x, y, returnType = TRUE)) ) convertIDs(x, y, ...)
	else x
}

#' Converting Gene or Probeset IDs
#'
#' The \pkg{CellMix} package implements methods for converting markers identifiers 
#' cross-platform and cross-organism, to facilitate deconvolution analysis to be carried 
#' out using multiple independent sources of data, e.g., use a marker gene list obtained 
#' from data on one platform to deconvolve gene expression data generated on another 
#' platform.
#' 
#' The function \code{convertIDs} provides the main interface to convert genes/probeset ids 
#' into IDs compatible with another given data.
#' It is typically useful to convert built-in marker gene lists (see \code{link{cellMarkers}}).  
#' 
#' @details 
#' The identifier conversion functions and methods defined in the \pkg{CellMix} package 
#' can be seen as extending the existing framework defined in the \pkg{GSEABase} package,  
#' with the generic \code{\link[GSEABase]{mapIdentifiers}}.
#' 
#' @param object object whose identifiers are converted.
#' @param to specification of the type of identifiers to convert to.
#' @param from specification of the type of identifiers of \code{object}.
#' This is only neeeded when the source type cannot be inferred from \code{object} itself.
#' @param ... extra arguments to allow extension, which are passed down to the workhorse method
#' \code{convertIDs,character,GeneIdentifierType,GeneIdentifierType}.
#' See each method's description for more details.
#' 
#' @export
#' @inline
setGeneric('convertIDs', function(object, to, from, ...) standardGeneric('convertIDs'))
#' Apply the conversion to each element of the list.
setMethod('convertIDs', signature(object='list', to='GeneIdentifierType', from='GeneIdentifierType')
		, function(object, to, from, ...){
			cl <- match.call()
			sapply(object, function(x) eval(cl), simplify=FALSE)
		}
)

#' Applies \code{mapIdentifier} to each element in a list.
#' 
#' All arguments in \code{...} are passed to the subsequent calls to 
#' \code{\link[GSEABase]{mapIdentifiers}}.
#' 
#' @param what identifiers to map
#'
#' @rdname convertIDs
#' @export
#' @inline
setMethod('mapIdentifiers', signature(what='list', to='GeneIdentifierType', from='GeneIdentifierType'),
	function(what, to, from, ..., verbose = FALSE){
		sapply(what, mapIdentifiers, to=to, from=from, ..., verbose=verbose, simplify=FALSE)
	}
)

#' Convert IDs from a MarkerList object.
#' 
#' In this case, argument \code{unlist} indicates if the result should be a simple list 
#' containing the mapping (a list) for each cell type or a \code{\linkS4class{MarkerList}} 
#' object (default).
#' 
#' @param verbose a logical or integer that sets the vverbosity level.
#' @param nodups specifies if marker identifiers that are duplicated across 
#' cell types in the result should be removed (\code{TRUE}) or not (\code{FALSE}).
#' If \code{NULL}, then duplicates are removed only if there were no duplicates in 
#' the source object.
#' 
#' @examples
#' 
#' # load a marker list from the registry
#' m <- MarkerList('IRIS')
#' summary(m)
#' head(m[[1]])
#' 
#' # convert Entrez gene ids to Affy probeset ids chip hgu133b 
#' m2 <- convertIDs(m, 'hgu133b.db', verbose=2)
#' summary(m2)
#'  
setMethod('convertIDs', signature(object='MarkerList', to='GeneIdentifierType', from='GeneIdentifierType'), 
		function(object, to, from, verbose=FALSE, nodups=NULL, unlist=TRUE, ...){
			
			# set verbosity level
			if( !missing(verbose) ){
				ol <- lverbose(verbose)
				on.exit( lverbose(ol) )
			}
			verbose <- lverbose()
			
			# logging functions
			log <- getLogger()
			vmessage <- log$message
			lmessage <- log$lmessage
			#
			
			# check for duplicated markers present in the original list
			hadDups <- !anyDuplicated(object)
			askDups <- !is.null(nodups)
			nodups <- if( askDups ) nodups else hadDups
			mapRes <- NULL
			nINI <- nmark(object)
			
			vmessage("# Converting ", nmark(object), " markers ", .GeneIdType_asString(from, to), ' ... ', appendLF=verbose>2)
			mnames <- marknames(object)
			gmap <- convertIDs(mnames, to=to, from=from, unlist = FALSE, ..., verbose=verbose-2)
			if( !is.list(gmap) )
				stop("convertIDs - Unexpected error: mapping of ", length(mnames)," marker names [", str_out(mnames), "] did not returned a list.")			
			toType <- attr(gmap, 'to')
			
			if( verbose ){
				stats <- .mtype(mappedkeysOnly(gmap[names(gmap) %in% mnames]), nINI)
				vmessage('OK [', stats, ']')
			}
			
			vmessage("# Processing ", nmark(object), " markers ", .GeneIdType_asString(from, to), ' ... ', appendLF=verbose>1)
			mids <- lapply(seq_along(object), function(i, ...){
				
				# add indentation
				oi <- log$nindent(1L)
				on.exit( log$nindent(oi, add=FALSE) )
				
				x <- object[[i]]
				type <- str_out(names(object)[i])
				if( !nchar(type) ) type <- str_c(i, '/', length(object))
				lmessage(2, "** Processing ids for ", type, " ... ", appendLF=FALSE)
				mnames <- marknames(x)
#				res <- convertIDs(mnames, to=to, from=from, unlist = FALSE, ..., verbose=verbose-2)
				res <- gmap[names(gmap) %in% mnames]
				
				if( !unlist ){
					lmessage(2, 'OK [', length(mappedkeys(res)), '/', length(x), ' (', .mtype(res), ')]')
					if( length(res) ) mapRes <<- c(mapRes, res)
					return(res)
				}
				
				## post-processing
				res <- mappedkeysOnly(res)
				if( length(res) ){
					res <- unlist2(res)
					ntot <- length(res)
					if( anyDuplicated(res) ){
						lmessage(3, "\n # Removing duplicated id(s) ... ", appendLF=FALSE)
						dups <- duplicated(res)
						res <- res[!dups]
						lmessage(3, "OK [dropped ", sum(dups) , "/", ntot," id(s)]")
					}
					
					mapRes <<- c(mapRes, res)
				}
				stats <- .mtype(res, length(x))
				# check if one should keep values
				if( hasValues(x) ){
					res <- setNames(x[names(res)], res)
				}
				
				lmessage(2, 'OK [', stats, ']')
				res
			} , ...)
	
			mids <- setNames(mids, names(object))
			
			# if only the map is requested: return it as a normal list/vector
			if( !unlist ){
				if( verbose ){
					stats <- .mtype(mappedkeysOnly(mapRes), nINI)
					vmessage('OK [', stats, ']')
				}
				return(mids)
			}
			object@.Data <- mids
			
			## post-processing
			lmessage(2, "# Checking for duplicated marker(s) across cell-types ... ", appendLF=FALSE)
			wd <- hasDuplicated(object, which=TRUE)
			nDUP <- sum(sapply(wd, length))
			if( !nDUP ){
				lmessage(2, 'OK')
			}else if(!nodups) {
				lmessage(2, if(askDups) 'SKIP' else 'OK', ' [', nDUP, ']')
			}else{ # remove duplicated probes
				n1 <- nmark(object)
				wd <- lapply(wd, '-')
				object <- object[wd]
				lmessage(2, "OK [dropped ", nDUP, "/", n1, "]")
			}
			
			# drop empty element
			object <- drop(object)
			
			if( verbose ){
				mapRes <- mapRes[mapRes %in% marknames(object)]
				stats <- .mtype(mapRes, nINI)
				vmessage('OK [', stats, ']')
			}
				
			# update annotation if necessary and possible
			geneIdType(object) <- NULL
			geneIdType(object) <- toType

			# return converted object			
			object
		}
)

keymap <- function(chip, src, dest, keys=NULL, reduce=TRUE){
	
	# load annotation package if necessary
	chip <- biocann_object(chip)
	
	kt <- keytypes(chip)
	if( !src %in% kt )
		stop("convertIDS - Invalid source key ['", src,"']: should be one of ", str_out(kt, Inf), '.')
	if( !dest %in% kt )
		stop("convertIDS - Invalid destination key ['", dest,"']: should be one of ", str_c(kt, Inf), '.')
	map <- AnnotationDbi::select(chip, keys=keys, cols=c(src, dest), keytype=src)
	
	if( !reduce ) return(map)
	gr <- factor(map[[src]], levels=keys)
	# somehow the column names are not guaranteed to be named as in argument cols
	cn <- colnames(map)
	akey <-  
			if( dest %in% cn ) dest 
			else{
				o <- which(cn!=src)
				if( length(o) == 0L ) 
					stop("Too few columns returned: could not choose other than ", src)
				if( length(o) > 1L )
					stop("Too many columns returned: could not choose between ", str_out(cn[o], Inf)) 
				cn[o]
			}
	# format as a list by key
	res <- by(map, gr, function(x) unique(as.character(x[[akey]])))
	sapply(res, identity, simplify=FALSE)	
}

#' This is the workhorse method that is eventually called by all other \code{convertIDs} 
#' methods.
#' The actual conversions are perforemd by \code{\link{mapIDs}}, to which are passed
#' all arguments in \code{...}, in particular, arguments \code{verbose} and \code{method}.
#' 
#' @param method mapping method, passed to \code{\link{mapIDs}}, that indicates
#' how to carry the mapping between the original and final identifier type. 
#' @param unlist logical that indicates if the result should be flatten, 
#' i.e. turned into a vector rather than a list -- using \code{\link{unlist2}}.
#' In this case, the vector's name then correspond to the source identifiers.
#' 
#' @examples
#' 
#' #----------------------------------------------
#' # 1. Conversion from biological IDs
#' #----------------------------------------------
#' # For this kind of IDs, a source annotation package can often be inferred 
#' # from the ID type, using regular expression patterns (e.g. "^ENS[0-9]+$" 
#' # identifies Ensembl gene IDs)
#'  
#' ids <- c("Hs.1", "Hs.2", "Hs.3")
#' # get Entrez gene IDs (based on annotation from the org.Hs.gene.eg package) 
#' convertIDs(ids, 'ENTREZID', 'org.Hs.eg.db', verbose=TRUE)
#' 
#' # map to other IDs 
#' convertIDs(ids, 'REFSEQ')
#' convertIDs(ids, 'ENSEMBL')
#' # convert across ogranism
#' convertIDs(ids, 'rat2302.db')
#' # get Affy probeset IDs for chip hgu133a
#' affy <- convertIDs(ids, 'hgu133a.db')
#' 
#' # assume we have a vector of IDs, e.g. Entrez gene ids
#' id <- c("673", "725", "10115")
#' # get associated probesets on chip hgu133a
#' convertIDs(id, 'hgu133a.db')
#' # get all associated probesets on chip hgu133a
#' convertIDs(id, 'hgu133a.db', method='all')
#' # same as a vector with duplicated names
#' convertIDs(id, 'hgu133a.db', method='all', unlist=FALSE)
#' # specification using ProbeAnnDbBimap objects
#' library(hgu133b.db)
#' convertIDs(id, 'hgu133a.db', hgu133bENTREZID, verbose=2)
#' 
#' #----------------------------------------------
#' # 2. Conversion from probeset IDs
#' #----------------------------------------------
#' # For this kind of IDs, a source annotation package is required, because it 
#' # cannot be easily inferred from the ID type.
#' 
#' # get Affy probeset IDs for chip hgu133b from ids for hgu133b 
#' convertIDs(affy, 'hgu133a.db', 'hgu133b.db')
#' # across organism
#' convertIDs(affy, 'hgu133a.db', 'rat2302.db')
setMethod('convertIDs'
		, signature(object='character', to='GeneIdentifierType', from='GeneIdentifierType'),   
	function(object, to, from, method='auto', unlist = TRUE, ...){
		
		# normalize to ensure both have annotations
		type <- .mapIdentifiers_normalize(from, to)
		from <- type[[1]]; to <- type[[2]]
		# extract annotation packages
		src <- .split_anndb(getAnnotation(from))
		dest <- .split_anndb(getAnnotation(to))
		
		## Compute mappings filling result list in place
		ids <- object
		ids <- ids[!duplicated(ids)]
		res <- NAmap(ids)
		
		.local <- function(x, y, ...){
			
			# for annotation packages or idtypes: map to all
			if( is.character(y) && is.null(names(y))  ){
				for( i in seq_along(y) ){
					
					if( length(y) == 1L ){
						# recursively fill missing mappings
						missing.id <- is.na(res) 
						lookup_ids <- ids[missing.id]
					}else{ 
						# cumulate/merge mappings
						lookup_ids <- ids
					}
					# setup local src and dest annnotations
					src <- from; 
					if( !is.null(x) ) annotation(src) <- x
					dest <- to; annotation(dest) <- y[i]
					# map ids
					lmap <- mapIDs(lookup_ids, src, dest, method=method, ...)
					# fill in global result
					#IMPORTANT: this is very sensitive piece of code
					# => careful on the impact of mapping
					if( length(y) == 1L ){
						# fill in mapped missing mappings
						mids <- lookup_ids[lookup_ids %in% names(lmap)]
						res[mids] <<- lmap[mids]
					}else{
						# merge mapped keys with previous mappings
						mapply(function(val, k){
							if( !is_NA(val) ){
								if( is_NA(res[[k]]) ) res[[k]] <<- val
								else res[[k]] <<- c(res[[k]], val[!val %in% res[[k]]])
							}
						}, lmap, names(lmap))
					}
				}
			}
#				else{ # dest is a mapping
#					# fill missing mappings
#					missing.id <- is.na(res)
#					res[missing.id] <<- mapIDs(ids[missing.id], x, dest, ...)
#				}
		}
		if( is.list(src) || (is.character(src) && length(src) > 1) ){
			if( identical(src, dest) )
				mapply(.local, src, dest, MoreArgs=list(...))
			else sapply(src, .local, y=dest, ...)
		}else .local(src, y=dest, ...)
		##
		
		## Post-process result
		res <- lapply(res, function(x) if(is.null(x)) NA else x)
		
		# flatten into a vector
		if( unlist ) res <- unlist2(res)
		# if possible attach attribute about the source and destination type
		if( !is_GeneIdentifierMap(to) )
			res <- ExposeAttribute(res, from=from, to=to, .VALUE=TRUE, .MODE='r')
		res
})


.convertIDs_matrix <- function(keys, to, from, ..., unlist, rm.duplicates){
	
	pids <- convertIDs(keys, to, from, ..., unlist = FALSE)
	if( !unlist ) return(pids)
	
	pids <- unlist2(pids)
	pids <- mappedkeysOnly(pids)
	if( anyDuplicated(pids) ){
		dups <- which(duplicated(pids))
		isdup <- which(pids %in% pids[dups])
		if( isFALSE(rm.duplicates) ){
			stop("Conversion returned ", length(isdup), " duplicated id(s): "
					, str_out(unique(pids[dups])))
		} else if( is.null(rm.duplicates) ){
			warning("Removed ", length(isdup), " duplicated row(s) after mapping: "
					, str_out(unique(pids[dups])))
			
			rm.duplicates <- TRUE
		}
		if( isTRUE(rm.duplicates) ){
			# remove duplicates
			pids <- pids[-isdup]
		}
	}
	pids
}

#' Convert the row names of a matrix into other identifiers.
#' 
#' In this case, argument \code{unlist} indicates if the converted ids should be 
#' used to subset the original \code{matrix} object, or returned directly returned
#' as a list.
#' 
#' @param rm.duplicates logical that how duplicated should be treated. 
#' \code{rm.duplicates=FALSE} does not allow any duplicated match and throws an error 
#' if any is present. If \code{TRUE} or \code{NULL} duplicates only the first match is kept, 
#' but a warning is thrown only when \code{NULL}.
#'  
setMethod('convertIDs'
		, signature(object='matrix', to='GeneIdentifierType', from='GeneIdentifierType')  
		, function(object, to, from, ..., unlist = TRUE, rm.duplicates=NULL){
			
			pids <- .convertIDs_matrix(rownames(object), to, from, ...
									, unlist=unlist, rm.duplicates = rm.duplicates)
			if( !unlist ) return(pids)
			# subset
			object <- object[names(pids), , drop=FALSE]
			rownames(object) <- pids
			# return
			object
		}
)
#' Convert the feature names of an ExpressionSet into other identifiers.
#' 
#' In this case, argument \code{unlist} indicates if the converted ids should be 
#' used to subset the original \code{ExpressionSet} object, or returned directly returned
#' as a list.
setMethod('convertIDs'
		, signature(object='ExpressionSet', to='GeneIdentifierType', from='GeneIdentifierType')  
		, function(object, to, from, ..., unlist = TRUE, rm.duplicates = NULL){
			
			pids <- .convertIDs_matrix(featureNames(object), to, from, ...
					, unlist=unlist, rm.duplicates = rm.duplicates)
			if( !unlist ) return(pids)
			
			# subset
			object <- object[names(pids), , drop=FALSE]
			featureNames(object) <- pids
			# update annotation if possible
			annotation(object) <- annotation(to)
			# return
			object
		}
)

#' Convert identifiers, inferring the type of origin from the object itself, 
#' but keep the annotation specification embedded in \code{from}. 
setMethod('convertIDs', signature(from='NullIdentifier'),  
		function(object, to, from, ...){
		
			# infer type from object
			from <- asGeneIdentifierType(object, annotation=annotation(from))
			if( isNullIdentifier(from) ){
				stop("convertIDS - Could not identify source type of ID keys: ", str_out(object)
						, ".\n    Please specify a valid ID type in argument `from`.")
			}
			# call next method
			convertIDs(object, to, from, ...)
		}
)

#' Convert identifiers, inferring the type from the specifications 
#' in \code{to} and \code{from}, eg., \code{to='ENTREZID'}, or \code{'UNIGENE'}.
#' If not specified in either \code{to} or \code{from}, the annotation is taken from 
#' \code{object}.
#' If \code{from} is missing, the source type is infered from \code{object} itself.
#'  
setMethod('convertIDs', signature(object='ANY'),  
	function(object, to, from, ...){
		
		if( !missing(to) && is(to, 'GeneIdentifierType') 
			&& !missing(from) && is(from, 'GeneIdentifierType') ){
			stop("ConvertIDs - No method is defined for object of class ", class(object))
		}
				
		# infer type from object
		objfrom <- asGeneIdentifierType(object)
		
		if( missing(from) ){
			from <- objfrom
			if( isNullIdentifier(from) ){
				stop("convertIDS - Could not infer source identifier type from keys: ", str_out(object)
						, ".\n    Please provide it in argument `from`.")
			}
		}
		
		# convert as GeneIdentifierType
		to <- asGeneIdentifierType(to)
		from <- asGeneIdentifierType(from)
		
		# use inferred type if it gives the id type
		if( !isAnnotationIdentifier(objfrom) && isAnnotationIdentifier(from) ){
			from <- asGeneIdentifierType(objfrom, annotation=annotation(from))
		}
		
		# use object annotation if necessary and possible
		if( hasAnnotation(object) && !hasAnnotation(to) && !hasAnnotation(from) ){
			annotation(from) <- annotation(to) <- annotation(object)
		}
		
		# infer destination type from specification
		convertIDs(object, to, from, ...)
	}
)

#' Convert identifiers using a given map or list of maps.
setMethod('convertIDs', signature(to='list', from='missing'),
	function(object, to, from, ...){
		if( isS4(to) ){
			convertIDs(object, asGeneIdentifierType(to), asGeneIdentifierType(object), ...)
		}else{
			convertIDs(object, GeneIdentifierMap(to), GeneIdentifierMap(), ...)
		}
	}
)
