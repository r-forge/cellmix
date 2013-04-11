# Class MarkerList
# 
# Author: Renaud Gaujoux
# Creation: 04 Mar 2011
###############################################################################

#' @include MarkerList-class.R
NULL

as.MarkerList <- function(x, ...){
	if( !is(x, 'MarkerList') ){
		new('MarkerList', .Data=x, ...)
	}else
		x
}

#' \code{isMarkerList} tests if an object is a \code{\linkS4class{MarkerList}} object.
#' 
#' @param x generally a \code{MarkerList} object or any object depending on the 
#' method. See their respective description for details.
#' 
#' @rdname markers 
#' @export
isMarkerList <- function(x){
	is(x, 'MarkerList')
}

.split_anndb <- function(x){
	if( isString(x) ){
		if( !nchar(x) ) ''
		else str_trim(strsplit(x, ',')[[1]])
	}else x
}

.glue_anndb <- function(x) paste(x, collapse=', ')

isNullIdentifier <- function(x){
	if( !is(x, 'GeneIdentifierType') ) x <- geneIdType(x)
	is(x, 'NullIdentifier')
}
isAnnotationIdentifier <- function(x){
	if( !is(x, 'GeneIdentifierType') ) x <- geneIdType(x)
	is(x, 'AnnotationIdentifier')
}

#' Extract name of the annotation package from a marker list.
setMethod('annotation', 'MarkerList', 
	function(object){
		.split_anndb(annotation(geneIdType(object))) 
	}
)
#' Sets the name of the annotation package associated to a marker list.
setReplaceMethod('annotation', signature('MarkerList', value='character'), 
	function(object, value){
		
		value <- .split_anndb(value)
		old_ann <- annotation(object)
		
		# early exit if nothing changes
		if( identical(value, old_ann) ) return(object)
		
		# check consistency if more than one map is passed 
		if( length(value) > 1L ){
			t <- sapply(value, function(m) idtype(biocann_object(m)))
			if( length(w <- which(t!=t[1L])) ){
				stop("Invalid annotation maps: their primary identifier should be of the same type"
								, " [", t[1L], ' != ', str_out(t[w], Inf), ']')
			}
		}
		
		# deal with the simple cases
		if( identical(value, '') || identical(old_ann, '') ){
			if( isNullIdentifier(object) ) geneIdType(object) <- value
			else if( value == '' && isAnnotationIdentifier(object) ) geneIdType(object) <- NullIdentifier()  
			else annotation(object@geneIdType) <- .glue_anndb(value)
			return( object )
		}
		
		# the annotations must share the same primary key
		newT <- idtype(biocann_object(value[1L]))
		oldT <- idtype(biocann_object(old_ann[1L]))
		if( !identical(newT, oldT) )
			stop("Primary annotation identifier types do not match: new type [", newT, " (", value[1L], ")]",
				" must be the same as the current type [", oldT, " (", old_ann[1L], ")].")
		
		# set new annotation
		annotation(object@geneIdType) <- .glue_anndb(value)
		# update annotation in geneIdType
		object
	}
)
#' Reset the annotation package associated to a marker list 
setReplaceMethod('annotation', signature(object='MarkerList', value='NULL'), 
	function(object, value){
		geneIdType(object) <- NULL
		object
	}
)
#' Returns the type of gene identifier used in a marker list.
setMethod('geneIdType', 'MarkerList', 
	function(object){
		object@geneIdType
	}
)
#' Sets the type of gene identifier used in a marker list to \code{\link{NullIdentifier}()}.
#' @param verbose logical not used.
setReplaceMethod('geneIdType', signature(object='MarkerList', value='NULL'), 
	function(object, verbose=FALSE, value){
		object@geneIdType <- NullIdentifier()
		object
	}
)
#' Sets the type of gene identifier used in a marker list to a given value.
setReplaceMethod('geneIdType', signature(object='MarkerList', value='character'), 
	function(object, verbose=FALSE, value){
		
		# special case for annotation package names
		if( all(is.annpkg(value)) ){
			geneIdType(object) <- AnnotationIdentifier(.glue_anndb(value))
			return(object)
		}
		
		value <- toupper(value)
		if( identical(value, 'NULL') ) value <- ''
		new_type <- asGeneIdentifierType(value, annotation=.glue_anndb(annotation(object)))
		geneIdType(object) <- new_type
		object
	}
)
#' Sets the type of gene identifier associated with a marker list to a given \code{\linkS4class{GeneIdentifierType}} 
#' object.
setReplaceMethod('geneIdType', signature(object='MarkerList', value='GeneIdentifierType'), 
	function(object, verbose=FALSE, value){
		# do nothing if not changing type
		from <- geneIdType(object)
		if( identical(value, from) ) return(object)
#		else if( isNullIdentifier(from) || isNullIdentifier(value) ){
			object@geneIdType <- value
			return(object)
#		}
#		# TODO: check this step
#		# do map
#		convertIDs(object, to = value, from = from, verbose = verbose)
	}
)

setAs('MarkerList', 'GeneSet', 
	function(from){
		
		# extract common slots
		s <- slotNames(from)
		s <- s[s %in% slotNames('GeneSet')]
		sl <- sapply(s, slot, object=from, simplify=FALSE)
		# add extra slots
		sl$geneIdType <- geneIdType(from)
		# extract all gene Ids
		ids <- unique( as.character(marknames(from, unlist=TRUE)) )
		do.call('GeneSet', c(list(ids), sl))
		
	}
)

######################################
# GSEABase classes coercions
######################################
.asGeneSetList <- function(object){
	# create common GeneSet template
	tpl <- as(object[0L], 'GeneSet')
	mapply(function(m, nm){
		initialize(tpl, geneIds=m, setName=nm)
	}, geneIds(object), names(object), SIMPLIFY=FALSE)
}

setAs('MarkerList', 'MarkerSetCollection', 
	function(from){
		new('MarkerSetCollection', .asGeneSetList(from))
	}
)

setAs('MarkerList', 'GeneSetCollection', 
	function(from){
		new('GeneSetCollection', .asGeneSetList(from))
	}
)

.showMarkerList <- function(object) {
	
	# cell type names
	cat("Types: ",
			str_c(selectSome(names(object), maxToShow=4), collapse=', '),
			" (total: ", length(object), ")",
			"\n", sep='')
	
	# type of list
	cat('Mode: ', mltype(object), "\n", sep='')
	
	## from GSEABase::showGeneSet
	cat("setName:", object@setName, "\n")
	mn <- marknames(object)
	cat("geneIds:",
			paste(selectSome(mn, maxToShow=4),
					collapse=", "),
			paste("(total: ", length(mn), ")\n",
					sep=""),
			sep=" ")
	show(geneIdType(object))
	show(object@collectionType)
	##
	
	# add geneValues if any
	v <- geneValues(object, unlist=TRUE)
	cat("geneValues:", 
			if( is.null(v) ) NA
			else paste(selectSome(v, maxToShow=4), collapse=", "),
			"\n")	
}

#' @rdname MarkerList-data
setMethod("show", "MarkerList",
	function(object) {
		cat("<object of class: ", class(object), ">\n", sep='')
		.showMarkerList(object)
		cat("details: use 'details(object)'\n")
	}
)
#' Shows more details than the regular \code{show} method.
#' @rdname MarkerList-data
setMethod("details", "MarkerList",
	function(object) {
		cat("<object of class: ", class(object), ">\n", sep='')
		.showMarkerList(object)
		object <- as(object, 'GeneSet')
		## from GSEABase:details,GeneSet
		cat("setIdentifier: ", setIdentifier(object), "\n", sep="")
		cat("description: ", description(object), "\n",
				if(nchar(longDescription(object))!=0 &&
						longDescription(object) !=  description(object)) {
					"  (longDescription available)\n"
				},
				"organism: ", organism(object), "\n",
				"pubMedIds: ", pubMedIds(object), "\n",
				"urls: ", paste(selectSome(urls(object), maxToShow=3),
						collapse="\n      "), "\n",
				"contributor: ", contributor(object), "\n",
				"setVersion: ", as(setVersion(object), "character"), "\n",
				"creationDate: ", creationDate(object), "\n",
				sep="")
		##
	}
)

#' Extracting Annotation from Objects
#' 
#' \code{hasAnnotation} tells if an object has some -- non empty -- attached annotation.
#' 
#' @param object an object
#' @param ... extra parameters (currently not used) 
#' 
#' @rdname annotation-utils
#' @export
hasAnnotation <- function(object, ...){
	!is.null( getAnnotation(object, ...) )
}
#' \code{getAnnotation} try extracting embedded annotations.
#' By default, it returns \code{NULL} if the object contains no annotation.
#' 
#' @param null logical that indicates if an empty character string should 
#' be return as \code{NULL}.
#' 
#' @rdname annotation-utils
#' @export
getAnnotation <- function(object, ..., null=TRUE){
	ann <- if( hasMethod('annotation', class(object)) ) annotation(object)
	else attr(object, 'annotation')
	if( null && isString(ann) && nchar(ann) == 0L ) ann <- NULL
	else if( !null && is.null(ann) ) ann <- '' 
	ann
}

#' Subsetting Marker Lists
#' 
#' @description
#' \code{MarkerList} objects have convenient subset methods, 
#' that allow to subset the list not only along the sets (first 
#' dimension) but also along the markers, which proves to be very useful
#' to subset markers list based on another set of identifier, e.g., the genes 
#' present in some expression data, or in another marker list.
#' 
#' The method \code{`[`} works performs basic strict subsets using integer, 
#' character and logical subsetting vectors.
#' 
#' @details
#' 
#' For \code{`[`}, argument \code{i} must be a \strong{valid} indexes,
#' i.e. a character vector containing only strings that match exactly some 
#' names of \code{x}, or an integer with indexes in the range of the length 
#' of \code{x}, or a logical vector not longer than \code{x}.
#' 
#' The second argument \code{j} is allowed not to match exactly, and may be 
#' used to subset the list against another set of identifiers, or 
#' limit each set to the top first markers.  
#' 
#' @inheritParams base::[
#' @param drop logical that indicates if one should drop empty cell types 
#' before returning the result -- as a \code{MarkerList}. 
#' 
#' @rdname subset
#' @export
#' @examples
#' 
#' # load IRIS markers and Abbas signature data 
#' m <- MarkerList('IRIS')
#' m
#' # only keep the markers present in a dataset
#' data(Abbas)
#' m[, featureNames(Abbas)]
#' # NB: this could also be done with subset(m, Abbas)
#' 
#' # take top 10 (or smaller)
#' summary(m[,1:10])
#' 
#' # take markers with associated score >= 0.9
#' nmark(m[m >= 0.9])
#' 
setMethod('[', 'MarkerList',
	function (x, i, j, ..., drop = FALSE) 
	{
		# quick exit for empty objects
		if( length(x) == 0L ){
			return(x)
		}
		
		if( !missing(i) ){
			#print("subset on i")
			nam <- names(x)
			x@.Data <- 
			if( is.logical(i) ){
				if( length(i) > length(x) )
					stop("Top logical index too long [", length(i),"]: length should me at most ", length(x))
				x@.Data[i]
			}else if( is.numeric(i) ){
				if( any(bad <- abs(i) > length(x)) )
					stop("Top integer index out-of-bound in MarkerList [", str_out(i[bad], 5), ']')
				x@.Data[i]
			}
			else{
				idx <- match(i, names(x), nomatch=0L)
				if( any(bad <- idx==0L) )
					stop("Top character index out-of-bound in MarkerList [", str_out(i[bad], 5), ']')
				i <- idx
				x@.Data[i]
			}
			#IMPORTANT: ensure that the names are also reordered
			x <- setNames(x, nam[i])
		}
		
		if( nmark(x) == 0L ){
			return(x)
		}
		# subset on the markers
		if( !missing(j) ){
			#print("subset on j")
			x@.Data <- lapply(x@.Data, function(y){
				if( is.logical(j) ){
					y[rep(j, length.out=length(y))]
				}else if( is.numeric(j) ){
					y[j[abs(j) <= length(y)]]
				}else{
					idx <- match(j, marknames(y, simplify=FALSE), nomatch=0L)
					y[idx[idx!=0L]]
				}
			})
		}
		
		if( drop ) drop(x) # setNames(x@.Data, names(x))
		else x
		
	}
)
#' This method subsets/reorders each set of marker using a subset specification found in 
#' 
#' @param match.names logical that indicates one should try to match the names of \code{x} 
#' and \code{i} before subsetting/re-ordering.
#' 
#' @rdname subset
#' @export
setMethod('[', signature(x='MarkerList', i='list'), 
	function (x, i, j, ..., match.names = TRUE, drop = FALSE) 
	{
		# quick exit for empty objects
		if( length(x) == 0L ){
			return(x)
		}
		
		# limit to the same set of cell types
		if( match.names && !is.null(names(i)) && !is.null(names(x)) ){
			idx <- match(names(x), names(i), nomatch=0L)
			i <- i[idx[idx!=0L]]
			x <- x[idx!=0L]
			if( length(x) == 0L ){
				return(x)
			}
		}
    # subset each set separately
		x@.Data <- mapply(function(y, idx){
			m <- marknames(y, simplify=FALSE)
			imn <-  marknames(idx, simplify=FALSE)
			# check for list of numeric used as integers
			if( is.null(imn) && is.numeric(idx) ){
        # TODO: add warning here
				imn <- as.integer(idx)
			}
			#
			if( is.logical(idx) ){ 
				if( is.null(imn) ) y[idx] # direct subset
				else if( is.integer(m) ){
					# subsetting an integer marker list:
					# it has no names so cannot match them
					# so we directly subset them.
					y[idx]
				}else{
					# TRUE markers only
					ti <- match(imn[idx], m, nomatch=0L)
					y[ti[ti!=0L]]
				}
				
			}else if( is.integer(imn) ){ # direct subset			
				if( length(imn) ) y[imn] else y
			}else{
				# use markernames
				ti <- match(imn, m, nomatch=0L)
				y[ti[ti!=0L]]
			}
		}, x@.Data, i, SIMPLIFY=FALSE)

		# further subset each element globally
		if( !missing(j) ) x <- x[, j, drop=FALSE]

		if( drop ) drop(x) # setNames(x@.Data, names(x))
		else x
	}
)

#' This method is equivalent to \code{x[i, , ..., match.names=FALSE]}, 
#' i.e. that each element of the marker list \code{x} are subset/reordered as 
#' its corresponding one in \code{j}, following the lists' order -- with
#' no attempt to match their elements names before subsetting.
#'
#' @rdname subset
#' @export
setMethod('[', signature(x='MarkerList', i='missing', j='list'), 
	function (x, i, j, ..., drop = FALSE){
		x[j, , ..., match.names=FALSE, drop = drop]
	}
)

#' Converting Marker Lists into data.frames
#' 
#' The S3 method \code{stack} allow to stack a \code{MarkerList} object 
#' into a \code{data.frame}, which may be a more convenient format to handle, 
#' depending on the computation one wants to perform.
#' 
#' @param x a \code{\link{MarkerList}} object for \code{stack}, a \code{data.frame}
#' as returned by \code{stack} for \code{unstack}.
#' @param ... extra arguments not used
#' 
#' @rdname stack
#' @importFrom utils stack
#' @S3method stack MarkerList
#' 
#' @examples
#'  
#' stack( rMarkerList(3, 3) )
#' stack( rMarkerList(3, 3, names=TRUE) )
#' stack( rMarkerList(3, 3, values=TRUE) )
#' 
stack.MarkerList <- function(x, ...){
	
	if( length(x) == 0L ) return( data.frame(values=numeric(), ind=character(), names=character()) )
	 
	# ensure there are names
	x <- addNames(x)
	
	df <- 
	if( hasValues(x) ){
		df <- utils:::stack.default(x)
		df$names <- marknames(x)
		df
	}else{
		df <- utils:::stack.default(x)
		df$names <- df$values
		df$values <- NA
		if( !anyDuplicated(df$names) )
			rownames(df) <- df$names
		df
	}

	# add a S3 class
	class(df) <- c("stackedMarkerList", class(df))
	df
} 

#' The method \code{unstack} performs the reverse operation of \code{stack}.
#' 
#' @rdname stack
#' @importFrom utils unstack
#' @S3method unstack stackedMarkerList
#' 
#' @examples
#' 
#' m <- rMarkerList(3, 3, values=TRUE)
#' sm <- stack(m)
#' identical(unstack(sm), m)
#' 
unstack.stackedMarkerList <- function(x, ...){
	res <- utils:::unstack.default(x, names~ind)
	
	# add values if necessary
	if( !all(is.na(x$values)) ){
		res <- setNames(sapply(seq_along(res), function(i){
			m <- res[[i]]
			j <- match(m, x$names)		
			setNames(x$values[j][x$ind[j]==names(res)[i]], m)
		}, simplify=FALSE), names(res))
	}
	# Build a MarkerList object
	MarkerList(res)
}

setAs('MarkerList', 'factor', function(from){
	
	nm <- marknames(from)
	factor(setNames(names(nm), nm), levels=unique(names(nm)))
	
})


.summaryMarkerList <- function(object){
	
	cat("<object of class MarkerList>\n")
	bkd <- nmark(object, TRUE)
	ntot <- if( length(object) == 0 ) 0 else sum(bkd)
	
	cat("Types:", length(object),
			if( length(object) > 0 ) 
				str_c("[", str_out(names(object), 3) ,"]"),
			"\n")
	# type of list
	cat('Mode: ', mltype(object), "\n", sep='')
	
	cat("Markers:", ntot, "\n")
	# show type of id
	idt <- idtype(object)
	cat("IDtype:", if(idt!='') idt else 'unknown'
			, if( ntot > 0 ) 
				str_c("[", str_out(marknames(object), 3, quote=idt!='_INDEX_') ,"]")
			, "\n")
	# show values
	if( hasValues(object) ){
		wv <- which(bkd!=0)
		cat("Values: [", str_out(object[[wv[1]]], 3), "]\n", sep='')
	}
	
	# show source annotation
	ann <- annotation(object)
	if( ann[1] != '' )
		cat("Source:",  paste(ann, collapse=", "), "\n")
	
	# return breakdown
	invisible(bkd)
}

#' @S3method print MarkerList
print.MarkerList <- function(x, ...) print.listof(x, ...)
			

setGeneric('summary', package='base')

#' Summary method for \code{MarkerList} objects
#' 
#' @param object a list of markers, usually a \code{MarkerList} object, or any 
#' object with an attached list of markers (cf. \code{attachMarkers} and 
#' \code{getMarkers}), depending on the method. 
#' See their respective description for details.
#' 
#' @return \item{summary}{ returns a numeric vector containing the number of markers for each cell-type}
#' @aliases summary,MarkerList-method
#' 
#' @usage \S4method{summary}{MarkerList}(object)
#' 
#' @rdname markers
#' @export
setMethod('summary', 'MarkerList', 
	function(object){
		
		res <- nmark(object, each=TRUE)
		attr(res, 'summary') <- capture.output(.summaryMarkerList(object))
		if( mltype(object, 'logical') ){
			attr(res, 'pass') <- nmark(object[object], each=TRUE)
		}
		class(res) <- c('MarkerList_summary', class(res))
		res
	}
)
#' @S3method print MarkerList_summary 
print.MarkerList_summary <- function(x, true.only=TRUE, ...){
	cat(attr(x, 'summary'), sep="\n")
	# show breakdowns
	pass <- attr(x, 'pass')
	if( !is.null(pass) ){
		cat("Breakdown (showing selected only):\n")
		show(head(pass, Inf))
	}
	if( sum(x) && (is.null(pass) || !true.only) ){
		cat("Breakdown:\n")
		show(head(x, Inf))
	}
}


## Annotate the genes specific to each cluster.
##
## This function uses the \code{annaffy} package to generate an HTML table from the probe identifiers. 
#setGeneric('annotate', function(x, annotation, ...) standardGeneric('annotate'))
#setMethod('annotate', signature(x='MarkerList', annotation='character'), 
#		function(x, annotation, filename='markers', outdir='.', name='Marker genes', ...)
#		{		
#			library(annaffy)
#			anncols<-aaf.handler()[c(1:3, 6:13)]			
#			
#			# add html suffix to filename if necessary
#			if( length(grep("\\.html$", filename)) == 0 ) filename <- paste(filename, 'html', sep='.')
#			
#			# for each cluster annotate the genes set		
#			print(summary(x))
#			x <- MarkerList(x)
#			lapply(names(x), function(n){
#						g <- x[[n]]
#						print(head(g))
#						if( length(g) == 0 ) return()
#						g <- as.character(g)
#						anntable <- aafTableAnn(g, annotation, anncols)
#						# generate HTML output
#						saveHTML(anntable, file.path(outdir,paste('M',n,filename)), title=paste(name, '[top', nrow(anntable),']'))
#					})
#			
#			# return nothing
#			invisible()
#		}
#)
#
#setMethod('annotate', signature(x='NMF', annotation='character'), 
#		function(x, annotation, ...)
#		{
#			s <- extractFeatures(x)
#			class <- .predict.nmf(t(s))
#			annotate(class, annotation=annotation, ...)
#		}
#)

#' Generating Marker Lists
#'
#' @description 
#' The functions documented here serve to generate marker lists and \code{\linkS4class{MarkerList}}
#' objects.
#' They are useful to develop and test deconvolution methods that use marker gene lists. 
#' 
#' \code{gmarkers} is the work horse generator function. 
#' It can generate "fixed" or random plain list of marker, optionnally with names and/or values.
#' 
#' @param n number of marker sets to generate, e.g., the number of cell types in the data.
#' @param k numeric (integer) value that specifies the number of markers per set.
#' If a negative value is specified, then it is interpreted as a growing factor for  
#' successive sets of markers, e.g., the default \code{k=-3} generates a first 
#' set of 3 markers (=3 x 1), a second set with 6 (= 2 x 3), a third with 9 (= 3 x 3), 
#' and so long.
#' @param names specification for names.
#' If a logical it indicates if names should be generated or not.
#' If a single string, it indicates the prefix to use for the names, the default 
#' being \dQuote{Marker.}.
#' If a character vector, then it specifies the names of each marker, and must then contain 
#' at least as many unique elements as markers in the whole generated list.
#' @param values specification for values.
#' If a logical, it indicates if a value should be generated and associated to each marker.
#' If a numeric vector, then it specifies the values to associate with each marker.
#' The names of this vector are currently not used.
#' If not shuffling (see argument \code{shuffle}) and \code{values} is a numeric vector, 
#' then it must contain at least as many unique elements as markers in the whole 
#' generated list.
#' @param shuffle logical that indicate if the markers should be shuffled randomly.
#' This is used only if fixed names or values have been specified.
#' 
#' @export
gmarkers <- function(n, k=-3, names=FALSE, values=FALSE, shuffle=FALSE){
	
	if( is.character(n) ){
		if( anyDuplicated(n) )
			warning("Duplicated types in argument `n`: reducing to unique types.")
		types <- unique(n)
		n <- length(types)
	} else if( is.numeric(n) ){
		types <- paste('Type', 1:n,sep='')
	} else
		stop("Invalid value for argument `n`: must be a single numeric or a character vector.")
	
	# convert from known potential types
	if( !isTRUE(names) && !isFALSE(names) ){
		cnames <- 
		if( is(names, 'ExpressionSet') ) featureNames(names)
		else if( !is.null(dim(names)) ) rownames(names)
		else names
		if( !is.character(cnames) )
			stop("Invalid value for argument `names` [",class(names),"]: could not extract character vector for names ")
		
		names <- cnames		
	}
	
	# min length for names and values
	minl <- 
	if( k >= 0 ) n * k 
	else (n+1) * n/2 * abs(k)
	
	prefix <- "Marker."
	if( is.character(names) ){
		
		if( isString(names) ){
			prefix <- names
			names <- TRUE
		}else{
			names <- names[!duplicated(names)]
			if( length(names) < minl )
				stop("Invalid value for argument `names`: character vector must contain at least ", minl, " unique elements.")
			# shuffle names
			if( shuffle ) names <- sample(names)
		}
	}
	if( is.numeric(values) ){
		if( !shuffle && length(values) < minl ){
			stop("Invalid value for argument `values`: numeric vector must contain at least ", minl, " values.")
		}
		# shuffle names
		if( shuffle && !is.character(names) ){
			values <- sample(values, minl, replace=length(values) < minl)
		}
	}
	
	res <- mapply(function(i, l){
				
		if( k >= 0 ){
			nm <- k
			st <- (i-1)*nm + 1
			lidx <- st:(st+nm-1)
		}else{
			k <- abs(k)
			nm <- i * k
			st <- (i-1)*i/2*k + 1
			lidx <- st:(st+i*k-1)
		} 
				
		names <- if( is.character(names) ) names[lidx]
		else if( isTRUE(names) || !isFALSE(values) ){
			paste(prefix, lidx, sep='')
		}
				
		if( isTRUE(values) ) setNames(runif(nm), names)
		else if( is.numeric(values) ) setNames(values[lidx], names)
		else if( is.null(names) ) rep(l, nm)
		else setNames(names, rep(l, nm))
			
	}, seq(n), types, SIMPLIFY=FALSE)	
	
	if( isFALSE(values) ){
		res <- unlist(res)
		if( shuffle ) res <- sample(res)
		res
	}
	else setNames(res, types) 
}


#' \code{gMarkerList} generates \code{MarkerList} objects.
#' 
#' This is a shortcut for \code{MarkerList(gmarkers(...))}
#' 
#' @param ... arguments eventually passed to \code{gmarkers}. 
#' 
#' @export
#' @rdname gmarkers
gMarkerList <- function(...) MarkerList(gmarkers(...))

#' \code{rmarkers} Generates a random plain list of markers.
#' 
#' It is a shorcut for \code{gmarkers(..., shuffle=TRUE)}.
#' @export
#' @rdname gmarkers
rmarkers <- function(...){
	gmarkers(..., shuffle=TRUE)	
}
#' Generates a random MarkerList object. 
#' 
#' This is a shortcut for \code{MarkerList(rmarkers(...))}
#' 
#' @export
#' @rdname gmarkers
rMarkerList <- function(...) MarkerList(rmarkers(...))


#' Flattening Marker Lists
#' 
#' \code{flatten} is an S4 generic function that flattens objects by unfolding 
#' their inner levels, like what \code{\link{unlist}} does for 
#' \code{\link{list}} objects.
#' 
#' @param object the object to be flattened
#' @param ... extra arguments to allow extension.
#' See each method's description for more details.
#' 
#' @seealso \code{\link{unlist}}
#' @inline
#' @export
setGeneric('flatten', function(object, ...) standardGeneric('flatten'))
#' For \code{\link{list}} objects this is equivalent to \code{\link{unlist2}}. 
#' 
#' @param use.names logical that indicates if the top level names, i.e. the 
#' names of list (e.g., the cell types) should be used instead of the 
#' names of each vector element.
#' 
#' @examples
#' 
#' l <- list(1,2,3,c(4,5,6))
#' flatten(l)
#' nl <- setNames(l, letters[1:length(l)])
#' flatten(nl)
#' flatten(nl, use.names=TRUE)
#' 
setMethod('flatten', 'list', function(object, use.names=FALSE){
		# unlist values only
		res <- unlist(object, use.names=FALSE)
		# extract exact names
		allnull <- TRUE
		nam <- lapply(seq_along(object), function(i){
			x <- object[[i]]
			if( use.names ){
				allnull <<- FALSE
				rep(names(object)[i], length(x))				
			}else if( is.null(names(x)) ) rep("", length(x))
			else{
				allnull <<- FALSE
				names(x)
			}
		})
		# set names to values
		setNames(res, if( !allnull ) unlist(nam, use.names=FALSE) else NULL)
	}
)
#' For \code{\link{MarkerList}} objects the conversion is similar to what 
#' \code{\link{unlist}} would do, but argument \code{use.names} is used 
#' slightly differently.
#' 
#' Secondly, the names of the main elements are append to the marker names 
#' (like \code{unlist} does if \code{use.names=TRUE}) only for marker lists that
#' contain numeric values (e.g. specificity scores). In the other cases, they 
#' are used either as values for character marker lists or as duplicated names 
#' for integer marker lists.
#' 
#' @param size an integer that specifies the desired size of the result vector, 
#' which is used only if the marker are identified by their index.
#' It is useful if one wants to use the result with data of known -- greater -- dimension.
#' 
#' @examples 
#' 
#' ml <- rMarkerList(3, names=TRUE)
#' flatten(ml)
#' flatten(ml, use.names=TRUE)
#' 
setMethod('flatten', 'MarkerList',
	function(object, use.names=FALSE, size=NA){
		
		if( length(object) == 0L ) return(list())
		
		# special case for integer lists
		if( is.integer(object[[1]]) ){
			object <- addNames(object, prefix=NULL)
			idx <- NULL
			mapply(function(l, x) idx <<- c(idx,setNames(x, rep(l, length(x)))), names(object), object)			
			res <- rep(NA, if( is_NA(size) ) max(idx) else size)
			res[idx] <- names(idx) 
			return(res)
		}
		
		use.names <- use.names & !is.null(names(object))
		res <- NULL
		sapply(seq_along(object), function(i){
			val <- object[[i]]
			if( use.names ){
				name <- names(object)[i]
				if( hasValues(val) ) names(val) <- paste(name, names(val), sep='.')
				else val <- setNames(rep(name, length(val)), val)
			}else if( !hasValues(val) )
				names(val) <- NULL
				
			res <<- c(res, val)
		})
		res
	}
)

#' @examples 
#' # get marker Ids
#' m <- rMarkerList(3, names=TRUE)
#' geneIds(m) 
setMethod('geneIds', 'MarkerList',
	function(object){
		marknames(setNames(object@.Data, names(object)), unlist=FALSE) 
	}
)
#' @export
setReplaceMethod('geneIds', signature(object='MarkerList', value='list'),
	function(object, value){
		object@.Data <- value
		object
	}
)
#' Returns the values embedded in a \code{MarkerList} object.
#' 
#' @param unlist a logical that indicates if the values should be
#' returned as a list (\code{FALSE}) or as a vector.
#' 
#' @examples
#' # get attached values
#' m <- rMarkerList(3, values=TRUE)
#' geneValues(m)
#'  
setMethod('geneValues', 'MarkerList',
	function(object, unlist=FALSE){
		l <- 
		if( hasValues(object) ) setNames(object@.Data, names(object))
		else return(NULL)

		if( unlist ) unlist2(l)
		else l
	}
)

#' \code{hasValues} tells if a marker list contains numeric values, e.g. specificity scores.
#' 
#' @rdname MarkerList-data
#' @export 
hasValues <- function(object, ...){
	
	if( is(object, 'GeneSet') ) is(object, 'GeneValueSet')
	else if( is.list(object) )
		length(object) != 0L && (storage.mode(object[[1]])=='double' ) #|| is.factor(object[[1]]) )
	else
		length(object) != 0L && (storage.mode(object)=='double' ) #|| is.factor(object) )
}


#' \code{marknames} is similar to \code{geneIds} but by default returns the all 
#' marker names unlisted, i.e. as a character vector.
#' 
#' @param simplify logical that is use when for logical marker lists only, and indicates
#' if all marker names should be returned (\code{simplify=FALSE}), 
#' or only the ones associated with a \code{TRUE} value.
#' 
#' @rdname MarkerList-data
#' @inline 
#' @export
setGeneric('marknames', function(object, ...) standardGeneric('marknames') )
#' @export
setMethod('marknames', 'vector', 
	function(object, simplify=TRUE){
		
		if( is.integer(object) ) object
		else if ( length(object) == 0L ) character()
		else if( is.character(object) ) setNames(object, NULL)
		else if( hasValues(object) ) names(object)
		else if( is.logical(object) ){
			if( simplify ) object <- object[object]
			names(object)
		}
		else
			stop('marknames - Could not extract names from object [',str_out(object, quote=FALSE, use.names=TRUE),']')
	}
)
#'
#' @examples
#' marknames( rMarkerList(3) )
#' marknames( rMarkerList(3, names=TRUE) )
#' marknames( rMarkerList(3, values=TRUE) )
#' 
setMethod('marknames', 'list', 
	function(object, unlist=TRUE, ...){
		
		if( !length(object) ){
			if( unlist ) return( character() )
			else return( list() )
		}
		
		# return the marker list
		if( unlist ){
			unlist2(lapply(object, marknames, ...))
		}else sapply(object, marknames, ..., simplify=FALSE)
	}
)

#' Utility Functions for MarkerList Objects 
#' 
#' \code{dropvalues} Drops the values associated with each marker (e.g. tissue specificity score), 
#' and returns an object of the same type as the input object.
#' It is a shorcut for \code{\link{marknames}(object, unlist=FALSE)}.
#' 
#' @param object an R object, typically a list that contains marker data 
#' structured in a similar way as in \code{MarkerList} objects.
#' @param ... extra argument to allow extension
#' @param int.ok logical that indicates if integer marker names are valid,
#' or if one should use their names instead.
#' 
#' @rdname utils-MarkerList
#' @export
dropvalues <- function(object, ..., int.ok=TRUE){
	if( int.ok || !mltype(object, 'integer') ) marknames(object, unlist=FALSE)
	else sapply(object, names, simplify=FALSE)
}

#' \code{mltype} returns the type of marker list or compare it with a given type:
#' \itemize{
#' \item character: elements are marker names/ids;
#' \item numeric: elements are numeric vectors (e.g., scores), named with the marker names;
#' \item index: elements are integer vectors corresponding to indexes in some reference
#' data;
#' \item logical: elements are logical vectors, named with marker names, which
#' result from logical operations on a \code{MarkerList} object.
#' }
#' 
#' @param type type to compare with.
#' 
#' @rdname utils-MarkerList
#' @export
mltype <- function(object, type=c('character', 'numeric', 'integer', 'logical')){
	
	# define type
	mlt <- if( !length(object) ) NA
	else class(object[[1L]])

	# return type or compare with reference
	if( missing(type) || is.null(type) ) mlt
	else{
		type <- match.arg(type)
		identical(mlt, type)
	}
}


setGeneric('.atrack', package='NMF')
#' Heatmap Annotation Track for MarkerList Objects
#' 
#' This method allows to add annotation tracks in heatmaps produced by
#' \code{\link{aheatmap}}, to highlight the position of markers, when 
#' plotting either the global expression values or cell type-specific 
#' signatures, whether measured or estimated.
#' 
#' @param object a MarkerList object
#' @param data the reference data to annotate
#' @param view a character string that specifies how markers
#' should be annotated:
#' \describe{
#' \item{\dQuote{single}}{markers are shown in a single track}
#' \item{\dQuote{split}}{each cell type is shown in a separate track}
#' \item{\dQuote{predict}}{one track for each column in \code{data} is added, 
#' and markers are shown on the track associated with the column corresponding to 
#' its maximum value.}
#' }
#' @param ... extra parameters passed to .atrack
#'  
#' @export
#' 
#' @examples
#' 
#' # load IRIS markers and the Abbas signature matrix
#' m <- MarkerList('IRIS')
#' data(Abbas)
#' aheatmap(Abbas, annRow=m)
#' 
#' 
setMethod('.atrack', 'MarkerList', 
	function(object, data=NULL, ..., view=c('single', 'split', 'predict')){
		
		# Currently the values are not supported
		object <- dropvalues(object)
		
		view <- match.arg(view)
		ann <- 
		switch(view
			, split = {
				ann <- sapply(seq_along(object), function(i){
					m <- object[i]
					nm <- marknames(m)
					factor(setNames(rep(names(object)[i], length(nm)), nm))
					#factor(flatten(m, use.names=TRUE))
				}, simplify=FALSE)
				setNames(ann, rep('Markers', length(ann)))
			}
			, predict = {
				
				refdata <-
				if( is.nmf(data) ) basis(data)
				else if( is(data, 'ExpressionSet') ) exprs(data)
				else data
		
				p_ma <- as.numeric(as.character(NMF:::.predict.nmf(refdata)))
				mi <- matchIndex(object, refdata)
				# build a vector that flags the marker rows with there type
				ma.type <- rep(NA, nrow(refdata))
				sapply(names(mi), function(n) ma.type[mi[[n]]] <<- n)
				# build separate annotation vector for each component
				rann <- sapply(1:ncol(refdata), function(i){
							x <- setNames(rep(NA,nrow(refdata)), rownames(data));
							# flag the markers predicted to be on the i-th component
							x[p_ma==i] <- ma.type[p_ma==i];
							#factor(x, levels=names(mi))
							x
#				  x1 <- x
#				  x1[ x==names(sm)[i] ] <- 'Agree'
#				  x1[ x!=names(sm)[i] ] <- 'Disagree'
#				  x1
#				  x[ x==names(sm)[i] ] <- paste('_',names(sm)[i], sep='')
				})
				rann <- as.data.frame(rann)
				#str(rann)
				rownames(rann) <- rownames(data)
				colnames(rann) <- rep("Markers", ncol(rann))
				rann
			}
			, single = list(Markers=factor(flatten(object, use.names=TRUE)))
		)
		.atrack(ann, data=data, ...)
	}
)

#' Creating Mapping from Marker Lists
#' 
#' @description
#' \code{matchIndex} match indexes contained in an object against some other reference data,  
#' returning an object that is interpretable relatively to the reference data, 
#' e.g. retreive the row index of a data matrix contained in a marker list.
#' 
#' Methods are provided for \code{MarkerList} objects to match against 
#' character vectors, \code{matrix}, \code{\link[Biobase]{ExpressionSet-class}} 
#' or \code{NMF} objects, etc..
#' 
#' @param object to match against the reference data
#' @param data reference data against which to match the object, 
#' i.e. retrieve the indexes from.
#' @param ... extra argument to allow extension   
#' 
#' @inline
#' @export
setGeneric('matchIndex', function(object, data, ...) standardGeneric('matchIndex'))
#' This is the workhorse method for \code{MarkerList} objects, that is eventually called 
#' by all other methods.
#' It matches marker names against names in a character vector and returns the corresponding 
#' sets of indexes for each marker set (cell type). 
#' 
#' @param unlist a logical that indicates if the indexes should be returned as a vector 
#' (\code{TRUE}) or another \code{MarkerList} object.
#' @param no.match a logical that specifies if the result should keep the 
#' markers that cannot be found in the reference object.
#' If \code{TRUE}, 0L values are used for unmatched markers.
#' Default is \code{FALSE}.
#' @param keep.names a logical that indicates if the markers' names should be kept
#' or removed (\code{FALSE}). In this case, the names of the sets of each markers are 
#' used -- and repeated.
#' @param verbose verbosity level, usually \code{TRUE} or \code{FALSE}.
#' 
#' @export
#' 
#' @examples
#' 
#' # load IRIS markers and Abbas signatures
#' m <- cellMarkers('IRIS')
#' data(Abbas)
#' summary(m)
#' ma <- matchIndex(m, Abbas)
#' summary(ma)
#' 
#' # get as a mapping vector
#' head(matchIndex(m, Abbas, unlist=TRUE))
#' # only set names
#' head(matchIndex(m, Abbas, unlist=TRUE, keep.names=FALSE))
#' 
#' # if keeping all one get 0-value indexes
#' summary(matchIndex(m, Abbas, no.match=TRUE))
#' 
setMethod('matchIndex', signature(object='list', data='character'),
	function(object, data, unlist=FALSE, no.match=FALSE, keep.names=TRUE, verbose=FALSE){
		
		# extract marker names from the list 
		lm <- dropvalues(object)
		# convert id if necessary
#		if( !missing(annotation) )
#			lm <- convertIDs(object, annotation)
		
		if( verbose ){
			message("Matching ", mltype(lm), " marker list against ", length(data), " strings [", str_out(data), ']'
					, " ... ", appendLF=FALSE)
		}
		# convert character index into integer index
		if( all(sapply(lm, is.character)) ){
			
			lm <- sapply(lm, function(x){
				idx <- match(x, data, 0L)
				if( keep.names ) idx <- setNames(idx, x)
				if( !no.match ) # remove unmatched markers
					idx[idx!=0L]
				else
					idx
			}, simplify=FALSE)
		}else if( mltype(lm, 'integer') ){
			lm <- sapply(lm, function(x) setNames(x, data[x]), simplify=FALSE)
		}
		
		if( verbose ){
			message('OK [', nmark(lm), '/', nmark(object), ' matche(s)]')
		}
		
		if( unlist ) unlist2(lm)
		else lm
	}
)
#' Match markers against the row names of a data matrix.
#' @export
setMethod('matchIndex', signature(object='list', data='matrix'),
	function(object, data, ...){
		
		if( is.null(rownames(data)) ){
			if( !all(sapply(object, is.integer)) )
				stop("matchIndex - data matrix has no rownames: cannot convert the marker list into an index list.")
			
			object
		}else{
			matchIndex(object, rownames(data), ...)
		}
	}
)
#' Match markers against the basis matrix of an \code{NMF} model object.
#' 
#' See \code{matchIndex,MarkerList,matrix-method}.
#' @export
setMethod('matchIndex', signature(object='list', data='NMF'),
	function(object, data, ...){		
		matchIndex(object, basis(data), ...)
	}
)
#' Match markers against the expression matrix of an \code{ExpressionSet} data object.
#' 
#' See \code{matchIndex,MarkerList,matrix-method}.
#' @export
setMethod('matchIndex', signature(object='list', data='ExpressionSet'),
	function(object, data, ...){
		matchIndex(object, exprs(data)
#			, annotation = if( missing(annotation) ) annotation(data) else annotation
			, ...)
	}
)
#' Match markers against the keys of a chip annotation database.
#' @export
setMethod('matchIndex', signature(object='list', data='ChipDb'),
	function(object, data, ...){
		pkg <- data$packageName
		matchIndex(object, keys(data), ...)
#		if( hasArg(annotation) ) matchIndex(object, keys(data), ...)
#		else matchIndex(object, keys(data)
#				, annotation = if( !is.null(pkg) && pkg != '' ) pkg else data, ...)
	}
)
#' Match markers against the keys of a probe annotation database.
#' @export
setMethod('matchIndex', signature(object='list', data='ProbeAnnDbBimap'),
	function(object, data, ...){		
		matchIndex(object, keys(data), annotation=data, ...)
	}
)
## #' The default method tries to extract marker data from \code{object} 
## #' using \code{\link{getMarkers}} and re-call \code{matchIndex} on it.
## #'   
## #' @export
# setMethod('matchIndex', signature(object='list', data='ANY'),
# 	function(object, data, ...){
# 		if( is.null(ma <- getMarkers(object)) ) return()
# 		matchIndex(ma, data, ...)
# 	}
# )
## #' @export
# setMethod('matchIndex', signature(object='list', data='missing'),
# 	function(object, data, ...){
# 		matchIndex(object, object, ...)
# 	}
# )
setMethod('matchIndex', signature(object='ANY', data='ANY'),
	function(object, data, ...){
		.Defunct(msg='matchIndex(ANY, ANY) - Not implemented')
	}
)
setMethod('matchIndex', signature(object='ANY', data='missing'),
	function(object, data, ...){
		.Defunct(msg='matchIndex(ANY, missing) - Not implemented')
	}
)

#' \code{attachMarkers} attaches a marker list to an \code{NMF} object, by storing it in an 
#' attribute named \code{'markers'}. It actually stores the results from 
#' \code{MarkerList(x)}, that is a \code{MarkerList} object.
#' 
#' @return \item{\code{attachMarkers}}{returns \code{object} with an attribute 
#' named \code{'markers'} that contains the result of \code{MarkerList(x)}.}
#' 
#' @export
#' @rdname markers
attachMarkers <- function(x, object){
	
	if( !is.null(object) ){
		object <- MarkerList(object)
		n <- if( is.nmf(x) ) nbasis(x) else ncol(x)
		if( n != length(object) )
			stop("attachMarkers - Invalid marker list: lengths do not match [nbasis(", nbasis(x),") != length(", length(object), ")].")
	}
	
	attr(x, 'markers') <- object
	x
}

#' \code{getMarkers} extracts the marker list from the attribute named \code{'markers'}
#' that will generally contain a \code{MarkerList} object.  
#'  
#' @param err.notfind a character vector that will prefix an error message if no
#' \code{x} has no marker list attached. The error is only thrown if 
#' \code{err.notfind} is not \code{NULL}. 
#' @param unlist a single logical specifying if the result should be unlisted
#' Default is \code{FALSE}. 
#'  
#' @return \item{\code{getMarkers}}{returns the attached marker list}
#'  
#' @rdname markers
#' @export
getMarkers <- function(x, err.notfind=NULL, unlist=FALSE){
	if( isNMFfit(x, recursive=FALSE) )	x <- fit(x)
	# load attribute 'markers'
	if( is.null(res <- attr(x, 'markers')) && is.nmf(x) )
		res <- attr(basis(x), 'markers')
	
	# throw an error if necessary
	if( !is.null(err.notfind) && is.null(res) )
		stop(err.notfind, " - could not find the marker list.")
	if( unlist )
		res <- unlist2(res)
	
	# return res
	res
}

#' \code{has.markers} tells if an object has an attached marker list, i.e. an 
#' attribute named \code{'markers'}.
#' 
#' @return \item{\code{has.markers}}{returns a logical}
#' 
#' @rdname markers
has.markers <- function(x) !is.null( getMarkers(x) )

#embedMarkers <- function(model, markers){
#	
#	stopifnot( length(markers) == nbasis(model) )
#	
#	# ensure the consistency of the identification of basis 
#	if( !is.null(names(markers)) && is.null(basisnames(model)) ) # transfer names from markers to model
#		basisnames(model) <- names(markers)
#	else if( is.null(names(markers)) && !is.null(basisnames(model)) )# transfer names from model to markers
#		names(markers) <- basisnames(model)
#	else if( !is.null(names(markers)) && !is.null(basisnames(model)) ){
#		# check that the names are the same and unique and enforce order 
#		stopifnot( !anyDuplicated(names(markers)) )
#		stopifnot( !anyDuplicated(basisnames(model)) )
#		idx <- match(basisnames(model), names(markers), 0L)
#		stopifnot( all(idx>0) ) 
#		# ensure the model and markers are in the same order if names are present	
#		markers <- markers[idx]
#	}
#	
#	# attach the markers to the model
#	model <- attachMarkers(model, as.MarkerList(markers))
#	
#	# return the model
#	model
#}

#' Enforcing Marker Block Patterns
#' 
#' This method enforces marker block patterns on a numeric matrix, which typically
#' is the basis matrix of the NMF model being estimated.
#' 
#' The block patterns are defined by a list of markers whose expression is cell-type
#' specific.
#' The expression value \eqn{e_{ij}}{e_ij} for a given marker \eqn{i}, on its 
#' respective cell-type \eqn{j} is forced to a given value; while its expression 
#' \eqn{e_{ik}}{e_ik} on any other cell-type \eqn{k} is forced to 0 or to be a
#' certain times lower than \eqn{e_{ij}}{e_ij}:
#' \eqn{e_{ik} = \min(\frac{e_{ij}}{f} - \epsilon, e_{ik})}{e_ik = min(e_ij/f - EPS, e_ik)}.  
#'    
#' 
#' @param object a numeric matrix or an \code{\link[NMF]{NMF-class}} object
#' 
#' @param markers a numeric vector or the list of markers to enforce. 
#' This is usually an \code{MarkerList} object, but can be any list,
#' as long as it contains markers that can matched to rows in \code{object}.
#' If \code{markers} is a single numeric, then it specifies the number of markers
#' to enforce in each cell-type. The markers are enforced on the top -- basis -- 
#' rows of \code{object}. 
#' If \code{markers} is a numeric vector, then its length must be the number of 
#' columns (basis) in \code{object}.
#' 
#' @param ratio a single numeric value specifying the minimum expression 
#' fold change to enforce. If \code{NULL} or zero, then the expression values 
#' of each marker are forced to 0 on non-related cell types.
#' Default is \code{NULL}.
#' 
#' @param value a single numeric value specifying the expression value to enforce 
#' on the cell-type related to each marker. If \code{NULL} this value is left 
#' unchanged. 
#' Default is \code{NULL}. 
#' 
#' @param eps minimum value imposed for null entries.
#' @param attach logical that indicates if the markers should be attached
#' to the matrix data, as attibute \code{'markers'}.
#' 
#' @param ... extra arguments to allow extension, and passed down to the workhorse
#' method \code{enforceMarkers,matrix,list}.  
#' 
#' @return the input object with the marker block patterns enforced
#' @author Renaud Gaujoux
#' 
#' @inline
#' @export
setGeneric('enforceMarkers', function(object, markers, ...) standardGeneric('enforceMarkers'))
#' Enforce marker patterns specified as a list.
#' 
#' This is the workhorse method that is eventually called by all other methods. 
setMethod('enforceMarkers', signature(object='matrix', markers='list')
	, function(object, markers, ratio=NULL, value=NULL, eps=NULL
				, attach=FALSE){
		
		stopifnot( length(markers) == ncol(object) )
		
		# convert the marker into indexes
		markers <- matchIndex(markers, object)	
		idx.markers <- unlist(markers, use.names=FALSE)		
		stopifnot( !anyDuplicated(idx.markers) )				
		
		# reorder the markers in the same order as the columns if necessary
		if( !is.null(names(markers)) && !is.null(colnames(object)) ){
			stopifnot( setequal(names(markers), colnames(object)) )
			markers <- markers[colnames(object)]
		}else if( is.null(colnames(object)) ){ # set the colnames only if none were defined
			colnames(object) <- head(names(markers), ncol(object))
		}else if( is.null(names(markers)) ){
			names(markers) <- head(colnames(object), length(markers))
		}
		# apply the constraints to the matrix (make sure to work on a copy)
		object <- neq.constraints.inplace(object, markers, ratio=ratio, value=value, copy=TRUE)
		
		# use given eps if necessary
		if( !is.null(eps) ){
			object[object==0] <- eps
		}
			
		if( attach )
			object <- attachMarkers(object, markers)
		
		# return results
		object
	}
)
#' Enforce a given number of marker patterns  
setMethod('enforceMarkers', signature(object='matrix', markers='numeric')
		, function(object, markers, ...){
			
			if( length(markers) == 1 )
				markers <- rep(markers, ncol(object))
			
			if( length(markers) != ncol(object) )
				stop("enforceMarkers - invalid length for argument `markers`: should be equal to the number of columns of argument `object`")
			
			# define marker list from breakdown
			l <- unlist(mapply(rep, colnames(object) %||% 1:ncol(object), markers, SIMPLIFY=FALSE))
			markers <- split(seq_along(l), l)
						
			enforceMarkers(object, markers, ...)
		}
)
#' The method for an \code{NMF} object enforces the markers on its basis matrix 
#' and sets its basis names to the marker list's names, if none were already defined. 
setMethod('enforceMarkers', signature(object='NMF', markers='ANY')
		, function(object, markers, ...){
			# enforce the markers on the basis matrix
			w <- enforceMarkers(basis(object), markers, ...)
			# get attached markers if any
			m <- getMarkers(w)
			# remove attached markers from basis: put them on the model itself
			w <- attachMarkers(w, NULL)
			basis(object) <- w
			# attach markers to the model
			object <- attachMarkers(object, m)
			
			# update the basisnames, in case these were set by enforceMarkers
			basisnames(object) <- colnames(basis(object))
			
			# return enforced object
			object
		}
)

#' @S3method subset MarkerList
subset.MarkerList <- function(x, subset, select, filter=NULL, skip=0L
			, total = FALSE
			, invert=FALSE, fixed=TRUE, ignore.case=FALSE, verbose=FALSE, ...){
	
	orig <- x
	# FIRST: select only the required elements of the list
	if( !missing(select) ){
		if( is(select, 'NMF') )
			select <- basisnames(select)
		
		if( verbose ){
			message("Select marker sets in ", str_out(select), ' ... ', appendLF=FALSE)
		}
		# subset and reorder the list itself
		x <- 
		if( is.logical(select) ) x[if( invert ) !select else select]
		else if( isString(select) && !fixed ){ # use grep				
			x[grep(select, names(x), invert=invert, ignore.case=ignore.case)]
		}else if( is.character(select) ){
			x[if( invert ) !(names(x) %in% select) else select]
		}else{
			x[ if( invert ) -select else select]
		}
		if( verbose ) message("OK [", length(x), "]")
	}
	
	# SECOND: subset each element of the list
	if ( missing(filter) && !missing(subset) && isNumber(subset) ){ 
		# shortcut for filtering
		filter <- subset
	}else if ( missing(filter) && missing(subset) && isNumber(total) ){ 
		# shortcut for filtering on total
		filter <- total
		total <- TRUE
	} else if( !missing(subset) ){
		
		if( is.annpkg(subset) ){
			x <- convertIDs(x, subset, ..., verbose=verbose)	
		}else{
			# convert if necessary
			x <- tryConvertIDs(x, subset, ..., verbose=verbose)
			#
			
			if( isMarkerList(subset) ){ # only keep/exclude markers found in MarkerList object
				subset <- marknames(subset)
			}
			# get the indexes relative to the subset and keep the unfound ones
			# NB: the matched keys are NMF: rownames, ExpressionSet: featureNames, character: the elements
			idx <- matchIndex(x, subset, no.match=TRUE, verbose=verbose)
			if( invert ){
				idx <- sapply(idx, function(id){
							id[id==0L] <- -1L
							id[id>0L] <- 0L
							id
						}, simplify = FALSE)
			}
			x <- x[idx != 0L] #mapply(function(xi, id) xi[id!=0], x, idx, SIMPLIFY=FALSE)
		}
	}
	
	# THIRD: skip the number of top-markers required
	if( skip > 0L )
		x <- sapply(x, function(xi) tail(xi, length(xi)-skip), simplify=FALSE)
	
	# FOUR: filter the number of markers required
	if( !is.null(filter) ){
		if( !total ) x <- sapply(x, head, filter, simplify=FALSE)
		else{
			s <- flatten_seq(x, n=filter)
			x <- attr(s, 'markers')
		}
	}
			
	# return the subset list
	x
}

#' \code{subset} subset a \code{MarkerList} object keeping only the markers 
#' that are present with a given reference set, which can be a character vector,
#' the rownames of a matrix or an \code{NMF} object, or the feature names of 
#' an \code{ExpressionSet}. 
#' The markers are matched using the function \code{matchIndex}, and gene 
#' identifier conversion is performed if necessary.
#' 
#' @param subset the reference set in which the markers are looked-up. It can be
#' a character vector, a matrix, or \code{NMF} or \code{ExpressionSet} objects. 
#' @param select a character, integer or logical vector specifying which set/cell-type to keep.
#' If a single string and \code{fixed=FALSE}, it is used as a regular expression to match 
#' against the set names.
#' @param filter a single integer specifying the maximum number of markers per cell-type 
#' @param skip a single integer specifying the number of top markers to skip in 
#' each cell type
#' @param total logical that inidicates if argument \var{filter} specifies the total
#' number of markers or the maximum number of markers per cell type.
#' @param invert logical used when \code{select} that indicates that one selects 
#' the set (cell type) that do not satisfy the condition defined by argument \code{select}. 
#' @param fixed logical used only when \code{select} is a single string to specify 
#' that it should be interpreted as a regular expression to match against the 
#' set names.
#' @param ignore.case logical used only when \code{select} is a single string to specify 
#' that the match should \strong{not} be case-sensitive.
#' @param verbose verbosity level, usually \code{TRUE} or \code{FALSE}.
#' @param ... extra arguments passed to \code{\link{convertIDs}} if necessary.
#' 
#' @return an \code{MarkerList} object
#' 
#' @rdname subset
#' @export
setMethod('subset', 'MarkerList', subset.MarkerList)

#' @importMethodsFrom BiocGenerics sapply
setGeneric('sapply', package='BiocGenerics')

#' Aplying Functions Along MarkerList Objects
#' 
#' Applies a given function to a \code{\linkS4class{MarkerList}} object.
#' 
#' @inheritParams BiocGenerics::sapply 
#' @param ASNEW logical that indicates if the returned object should 
#' have all data other than the list identical to \code{X}.
#' @param ARG.NAMES logical that indicates if the names of each element should be passed to \code{FUN}
#' as a second argument. 
#' In this case, the function is applie using \code{\link{mapply}}, with two vectorised arguments: the list itself and its names. 
#' @param simplify logical that indicates if the result should be simplified
#' if possible.
#' If \code{TRUE}, then the result is simplified as in \code{\link[base]{sapply}}.
#' If \code{FALSE}, then the result is not simplified, and returned as a 
#' \code{\linkS4class{MarkerList}} object if possible, or as a list otherwise.
#' If \code{NA}, then the result is not simplified, and always returned as a list.
#'  
#' @export
setMethod('sapply', 'MarkerList',
	function(X, FUN, ..., ASNEW=FALSE, ARG.NAMES = FALSE, simplify=TRUE, USE.NAMES=TRUE){
		
		# early exit if empty object
		if( !length(X) ){
			return( if(simplify) list() else X )
		}
		
		# do not simplify if replace is TRUE
		.simplify <- simplify
		if( ASNEW ) .simplify <- FALSE
		else if( is_NA(simplify) ) .simplify <- FALSE 
		
    if( !ARG.NAMES ){
		  res <- sapply(setNames(X@.Data, names(X)), FUN, ..., simplify = .simplify, USE.NAMES = USE.NAMES)
    }else{
      res <- mapply(FUN, X@.Data, names(X), MoreArgs = list(...), SIMPLIFY = .simplify, USE.NAMES = USE.NAMES)
    }
		if( ASNEW ){
			return( MarkerList(res) )
		} else if( isFALSE(simplify) ){
			X@.Data <- res
			if( !is(try( validObject(X), silent=TRUE), 'try-error') )
				return(X)
		}
		res
	}
)

getCurrentGeneric <- function(envir=parent.frame()){
	f <- get(".Generic", envir = envir)
	get(f, envir=envir)
}


setGeneric('Compare', package='methods')
#' Operations on MarkerList Objects
#' 
#' \code{Compare} compares all embedded values with a given fixed value.
#' This is useful to filter markers based on their associated scores.
#' 
#' @inheritParams methods::Compare
#' 
#' @export
setMethod('Compare', signature(e1='MarkerList', e2='numeric'),
	function(e1, e2){
		if( !hasValues(e1) && !mltype(e1, 'integer') ){
			warning("Comparing MarkerList with no embedded values or index to a numeric: returned object unchanged.")
			return(e1)
		}
		
		cmp <- getCurrentGeneric()
		sapply(e1, cmp, e2, simplify=FALSE)
	}
)

#' \code{rmDuplicated} remove marker identifiers that are duplicated either within 
#' or between sets.
#' Arguments in \code{...} are not used.
#' 
#' @export
#' @rdname utils-MarkerList
rmDuplicated <- function(object, ...){
	mn <- marknames(object)
	if( anyDuplicated(mn) ){		
		object <- subset(object, subset=mn[duplicated(mn)], invert=TRUE)
	}
	object
}

#' \code{hasDuplicated} checks for duplicated identifiers across sets.
#' 
#' @param which logical that indicates if a list of duplicated identifiers 
#' should be returned.
#' If \code{which=TRUE} and no duplicates exist then \code{NULL} is returned.
#' 
#' @rdname utils-MarkerList
#' @export
hasDuplicated <- function(object, which=FALSE){	
	
	if( !which ){
		mn <- marknames(object)
		return( anyDuplicated(mn) )
	}
	res <- sapply(seq_along(object), function(i){
		dups <- marknames(object[[i]]) %in% marknames(object[-i])
		which(dups)
	})
	
	res
}

#test2 <- function(mart){
#	qu <- "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' uniqueRows = '1' count = '0' datasetConfigVersion = '0.6' requestid= 'biomaRt'> <Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'affy_hg_u133a'/><Attribute name = 'affy_hg_u133b'/><Attribute name = 'entrezgene'/><Filter name = 'affy_hg_u133b' value = '242904_x_at,232165_at,201369_s_at,233127_at,204621_s_at,202861_at,222862_s_at,235213_at,235739_at,231798_at,204470_at,210118_s_at,217767_at,221463_at,210511_s_at,205922_at,203591_s_at,223809_at' /></Dataset></Query>"
#	host <- if( missing(mart) ) 'http://web.cbio.uct.ac.za/~renaud/bimart.php' else biomaRt:::martHost(mart)
##	host <- 'http://www.cnn.com'
#	RCurl::postForm(host, query=qu, style='HTTPPOST'
#	#RCurl::getURL(host 
#	, .opts=list(verbose=TRUE
##			, proxy = "campusnet.uct.ac.za"
#			, proxy = "localhost"
##			, proxyauth = "ntlm"
##			, proxypassword = "tdbtlg34@"
#			, proxyport = 8080
##			, proxy.transfer.mode       
##			, proxytype = 'ntlm'                   
##			, proxyusername = 'gjxren001'
##			, proxyuserpwd = "wf.uct.ac.za\\GJXREN001:tdbtlg34@"
##			, proxyuserpwd = "gjxren001@wf.uct.ac.za:tdbtlg34@ "
#			))
#}
#
#test <- function(){
#	
##	qu <- "<?xml version='1.0' encoding='UTF-8'?><!DOCTYPE Query><Query  virtualSchemaName = 'default' uniqueRows = '1' count = '0' datasetConfigVersion = '0.6' requestid= 'biomaRt'> <Dataset name = 'hsapiens_gene_ensembl'><Attribute name = 'affy_hg_u133a'/><Attribute name = 'affy_hg_u133b'/><Attribute name = 'entrezgene'/><Filter name = 'affy_hg_u133b' value = '242904_x_at,232165_at,201369_s_at,233127_at,204621_s_at,202861_at,222862_s_at,235213_at,235739_at,231798_at,204470_at,210118_s_at,217767_at,221463_at,210511_s_at,205922_at,203591_s_at,223809_at' /></Dataset></Query>"
##	postForm(biomaRt:::martHost(mart), query=qu, style='post')
##	aff <- fData(Abbas)$AFFY[1:50]
#	aff <- cellMarkers('Abbas', id='AFFY')
#	i <- c('hgu133a.db', 'hgu133b.db')
##	i <- c('affy_hg_u133a')
#	o <- 'ENTREZ'
#	convertIDs(aff, i, o)	
#}

#' \code{reverse} reverse is an S4 method to reverse an object.
#' For \code{MarkerList} object, it reverse the order of markers within each set seprately, 
#' and optionally the order of the set in the list.
#' 
#' @inline
#' @rdname utils-MarkerList
#' @export
setGeneric('reverse', function(object, ...) standardGeneric('reverse'))
#'
#' @param all a logical that indicates if the list itself, i.e., the sets should be 
#' reversed as well.
#' 
#' @examples
#' m <- gMarkerList(3, 4)
#' geneIds(m)
#' geneIds( reverse(m) )
#' geneIds( reverse(m, all=TRUE) )
#' 
setMethod('reverse', 'MarkerList',
  function(object, all=FALSE){
	# reverse each element
	object@.Data <- lapply(object@.Data, rev)
	# reverse top list as well if requested
	if( all ) object <- object[length(object):1]
	
	object
  }
)

#' \code{nmark} returns the number of markers for each cell-type. It accepts any 
#' object as an input, but is generally called on an \code{MarkerList} object 
#' or an object with an attached \code{MarkerList} object (see methods 
#' \code{getMarkers} and \code{attachMarkers}).
#' 
#' @aliases nmark,MarkerList-method nmark,ANY-method 
#' 
#' @export
#' @inline
#' @rdname utils-MarkerList
setGeneric('nmark', function(object, ...) standardGeneric('nmark') )
#' The method for lists, including \code{MarkerList} objects returns the total 
#' number of markers, or optionaly its breakdown per set of markers.
#' 
#' @param each a logical indicating if the number of marker in each set is to 
#' be returned or only the total number (default).
#' 
setMethod('nmark', 'list'
	, function(object, each=FALSE){
		n <- sapply(object, length)		
		
		# return number list or total number
		if( !length(object) ) 0L
		else if( each ) n
		else sum(n) 
	}
)
#' The default method tries to extract embedded marker data from \code{object}
#' using \code{\link{getMarkers}}, and returns the result of \code{nmark} applied 
#' to these data if present or 0L otherwise.
#' 
setMethod('nmark', 'ANY'
	, function(object, ...){
		n <- if( has.markers(object) ){
			ma <- getMarkers(object)
			if( is(ma, 'MarkerList') )
				nmark(ma, ...)
			else
				0
		}
		else
			0
	}
)

#' Returns the character vector of all marker names in the list.
setMethod('featureNames', 'MarkerList', function(object) marknames(object) )

setGeneric('drop', package='base')

#' \code{drop} drops empty sets of markers from a \code{MarkerList} object, 
#' as well as markers with \code{FALSE} values if the object is a logical 
#' marker list.
#' It returns the reduced \code{MarkerList} object.
#' 
#' @export
#' @rdname utils-MarkerList
#' @examples
#' m <- MarkerList( list(1:4, 10:15, integer(), 100) )
#' geneIds(m)
#' geneIds(drop(m))
#' 
setMethod('drop', 'MarkerList', 
	function(x){
		x <- x[nmark(x, each=TRUE)>0L]
		# drop FALSE markers in logical lists
		if( mltype(x, 'logical') )
			x <- x[, featureNames(x)]
		x
	}
)


#' @export 
setGeneric('reorder', package='stats')
#' Reordering Marker Lists
#' 
#' The method \code{reorder} for \code{MarkerList} objects allows to reorder 
#' each set of markers according to some auxiliary reference data, 
#' e.g, their gene expression in different cell types.
#' 
#' When \code{data} is a matrix-like object, a score is computed for each marker in two steps: 
#' the first step computes scores for each marker within each group of columns defined by the 
#' levels of \code{by} (function \code{rowMins});
#' the second step aggregates these scores into a single value for each marker (function \code{fun}).
#' 
#' The markers are then ordered wihtin their respective set, according to the score and 
#' argument \code{decreasing}.
#' 
#' @param x a \code{MarkerList} object
#' @param data reference data as a numeric vector, a matrix, or an \code{ExpressionSet} 
#' object.
#' These data must have names or rownames against which markers can be match.
#' @param by optional \code{factor} only used when \code{data} is a matrix-like object, 
#' to 
#' @param fun function to perform the second step in the computation of scores (cf. argument \code{by}).
#' @param by.fun function to perform the first step in the computation of scores (cf. argument \code{by}).
#' @param decreasing logical that indicates how the markers should be ordered according to their score.
#' 
#' @return a \code{MarkerList} object that also contain the computed scores. 
#' 
#' @seealso \code{\link{reverse,MarkerList-method}}
#' 
#' @export
#' @rdname sort
setMethod('reorder', 'MarkerList'
	, function(x, data, by=NULL, fun=max, by.fun=rowMins, decreasing=!missing(data)){
		
		# plain reordering
		if( missing(data) ){
			x <- sort(x, decreasing=decreasing)
			return(x)
		}
		
		X <- data
		x <- subset(x, X)
		if( nmark(x) == 0 ){
			warning("None of the markers in `x` were found in `data`: returned an empty MarkerList object."
					, "\n  Check compatibility of key types.")
			return(x)
		}
		# use expression matrix for ExpressionSet objects
		if( is(X, 'ExpressionSet') ) X <- exprs(X)
		
		mn <- dropvalues(x)
		# for matrix object: compute single score by row
		scores <- 
		if( is.matrix(X) ){
			if( is.null(rownames(X)) )
				stop("reorder - Invalid value for reference `X`: matrix object should have rownames.")
			scores <- X[unlist(mn),,drop=FALSE]			
			if( !is.null(by) ){
				scores <- rowApplyBy(scores, by, by.fun)
			}
			apply(scores, 1L, fun)	
		}else if( is.numeric(X) ){
			if( is.null(names(X)) )
				stop("reorder - Invalid value for reference `X`: vector object should have names.")
			X[unlist(mn)]
		}
		else 
			stop('reorder - Invalid value for argument `X`: expect a matrix-like object or a numeric vector.')		
		
		if( is.null(names(scores)) )
			stop("reorder - Unexpected error: score vector should have names.")
		
		# reorder within each type
		x@.Data <- lapply(mn, function(mnames){
			sort(scores[mnames], decreasing=decreasing)
		})

		# return x 
		x
	}
)


#' \code{sort} sorts each element of a marker list.
#'  
#' @S3method sort MarkerList
#' @rdname sort
sort.MarkerList <- function(x, decreasing=FALSE, ...){
	
	if( hasValues(x) )
		sapply(x, function(x) x[order(x, decreasing=decreasing)], simplify=FALSE)
	else{
		sapply(x, function(x) x[order(marknames(x), decreasing=decreasing)], simplify=FALSE)
	}
}

numprefix <- function(x, n=length(x), prefix=''){
	pat <- str_c("%0", nchar(as.character(n))+1, "i")
	paste(prefix, sprintf(pat,1:length(x)), sep='')
}

flatten_seq <- function(x, decreasing=FALSE, n=NULL){
	
	nm <- sapply(x, length)
	nmax <- max(nm)
	idx <- mapply(function(a, x){
				setNames(x, paste(numprefix(x, nmax), a, sep='_'))
			}, numprefix(names(x), prefix='T'), x)
	idx <- unlist2(idx)
	if( !is_NA(decreasing) )
		idx <- idx[order(names(idx), decreasing = decreasing)]
	
	# add some attributes
	gr <- factor(sub("^.*_(T[0-9]+)$", "\\1", names(idx)))
	bd <- summary(gr, maxsum=Inf)
	# limit to a certain total number of elements
	if( !is.null(n) ){
		n <- as.numeric(n)
		idx <- head(idx, n)
		bd <- summary(head(gr, n), maxsum=Inf)
		i <- 0L
		attr(idx, 'markers') <- sapply(x, function(x){ 
					i <<- i + 1L
					head(x, bd[i]) 
				}, simplify=FALSE)
	}
	attr(idx, 'breakdown') <- bd
	# return flattened object
	idx
} 

#' @importFrom BiocGenerics combine
#' @export
setGeneric('combine', package='BiocGenerics')
#' Combine markers from multiple cell types of a MarkerList object, based on groups 
#' defined by a given factor.
setMethod('combine', signature(x='MarkerList', y='factor')
	, function(x, y, ...){
		
		if( length(y) != length(x) )
			stop("combine - Could not combine types in MarkerList object: argument `y` <factor> must be of the same length as `x` <MarkerList>.")
		
		idx <- split(seq_along(x), y)
		# combine the markers from each group of cell types
		x@.Data <- lapply(idx, function(i){
			unlist2(x@.Data[i])
		})
		# use levels as new cell type names
		names(x) <- levels(y)
		# return modified x
		x
	}
)

#' Calls an NMF algorithm using a MarkerList object, whose length defines the factorization 
#' rank.
#' The actual MarkerList object is passed down in argument \code{data}.
#' 
setMethod('nmf', signature(x='MatrixData', rank='MarkerList'), 
	function(x, rank, method, ...){
		if( missing(method) ) method <- NULL
		nmf(x, length(rank), method=method, ..., data=rank)
	}
) 
 