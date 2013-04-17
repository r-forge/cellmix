# Methods for class ExpressionMix
# 
# Author: Renaud Gaujoux
# Created: 24 Jul 2012
###############################################################################

#' @include ExpressionMix-class.R
#' @include MarkerList-class.R
NULL

isExpressionSet <- function(x){
	is(x, 'ExpressionSet')
}

isExpressionMix <- function(x){
	is(x, 'ExpressionMix')
}

#' @rdname ExpressionMix
#' @export
setMethod('show', 'ExpressionMix', 
	function(object){		
		callNextMethod()
		# add information about the mixture
		if( nbasis(object) > 0 ){
			cat("Composition: ", str_out(basisnames(object)), " (", nbasis(object), " total)\n", sep='')
		}else cat("Composition: none\n")
	}
)

#' Dimensions in ExpressionMix Objects
#' 
#' @description
#' Similarly to \code{\linkS4class{NMF}} objects, \code{ExpressionMix} objects
#' have a "third" dimension, the number of underlying cell types, inherited 
#' from the embedded \code{NMF} model.
#' 
#' \code{dim} returns the dimensions of an ExpressionMix object. 
#' It returns a 3-length integer vector, containing the number of features, 
#' samples and components respectively.
#' 
#' @param x an \code{ExpressionMix} object
#' 
#' @export
#' @rdname ExpressionMix-dims
setMethod('dim', 'ExpressionMix', 
	function(x){
		c(selectMethod('dim', 'ExpressionSet')(x), Components=nbasis(x))
	}
)

#' \code{dimnames} returns the dimension names of an ExpressionMix object.
#' It returns a list with 3 elements: the feature names, sample names and the 
#' constituent cell/tissue names (i.e. the basis names of the underlying NMF model).
#' 
#' @seealso \code{\link{featureNames}}, \code{\link{sampleNames}}, 
#' \code{\link{basisnames}} 
#' @export
#' @rdname ExpressionMix-dims
setMethod('dimnames', 'ExpressionMix', 
		function(x){
			list(featureNames(x), sampleNames(x), basisnames(x))
		})

#' \code{dimnames<-} Sets the dimension names of an ExpressionMix object.
#' The replacement value must be a list containing the feature names, the samples names
#' (optional), and the constituent cell/tissue names (optional), all as character 
#' vectors of correct dimensions (i.e. compatible with the number of features, 
#' samples and constituents).
#' 
#' @param value replacement value
#' 
#' @export
#' @rdname ExpressionMix-dims
setReplaceMethod('dimnames', 'ExpressionMix', 
		function(x, value){
			if( !is.list(value) )
				stop("ExpressionMix::dimnames - Invalid value: must be a list.")
			
			if( length(value) == 0 )
				value <- NULL
			else if( length(value) == 1 )
				value <- c(value, list(NULL, NULL))			
			else if( length(value) == 2 ) # if only the two first dimensions reset the third one
				value <- c(value, list(NULL))
			else if( length(value)!=3 ) # check length of value
				stop("ExpressionMix::dimnames - invalid argument 'value' [a 2 or 3-length list is expected]")
			
			# update NMF rownames and colnames if necessary		
			if( hasBasis(x) ) dimnames(basis(x)) <- value[c(1,3)]
			if( hasCoef(x) ) dimnames(coef(x)) <- value[c(3,2)]		
			# update the ExpressionSet dimnames
			fn <- selectMethod('featureNames<-', 'ExpressionSet')
			x <- fn(x, value=value[[1]])
			fn <- selectMethod('sampleNames<-', 'ExpressionSet')
			x <- fn(x, value=value[[2]])
			# set names on `pure` if necessary
			if( length(x@pure) > 0 ) 
				x@pure <- setNames(x@pure, sampleNames(x))
			
			# return updated x
			x
		})

#' Acessing Pure or Mixed Sample Data  
#' 
#' \code{npure} returns the number of pure samples.
#' 
#' @param object an \code{\linkS4class{ExpressionMix}} object.
#' 
#' @export
npure <- function(object){
	p <- attr(object, 'pure')
	if( length(p) == 0 ) 0L
	else sum(p)
}

#' \code{lpure} returns the pure samples as a logical vector.
#' 
#' @rdname npure
#' @export
lpure <- function(object){
	p <- attr(object, 'pure')	
	if( length(p) == 0L ) setNames(rep(FALSE, ncol(object)), colnames(object))
	else p	
}

#' \code{wpure} returns the indexes of the pure samples.
#' 
#' @rdname npure
#' @export
wpure <- function(object) which(lpure(object))

#' \code{pureSamples} returns the data for the pure samples only.
#' 
#' @param drop logical that indicates if empty levels should be removed
#' from phenotypic factor covariates of the returned \code{ExpressionMix}
#' object.
#' 
#' @rdname npure
#' @export
pureSamples <- function(object, drop=TRUE){
	# subset
	object <- object[, lpure(object)]
	# drop levels from the phenotypic annotation data.frame if requested
	if( drop ) pData(object) <- droplevels(pData(object))
	# return result
	object
}

#' \code{mixedSamples} returns the data for the mixed samples only.
#' 
#' @rdname npure
#' @export
mixedSamples <- function(object, drop=TRUE){
	# subset
	object <- object[, !lpure(object)]	
	# drop levels from the phenotypic annotation data.frame if requested
	if( drop ) pData(object) <- droplevels(pData(object))
	# return result
	object
}

#' Shows which data is available in an \code{ExpressionMix} object.
#' 
#' This method invisibly returns a 3-length logical that indicates if expression, 
#' signature and proportion data are available in \code{object}.
#' 
#' @param show logical that indicates if something is to be displayed on the 
#' standard output or not.  
#' 
setMethod('showData', 'ExpressionMix', 
		function(object, show=TRUE){
			res <- c(exprs = nrow(exprs(object))!=0,
					signatures = hasBasis(object),
					proportions = hasCoef(object))
			if( show ){
				d <- dim(object)
				if( res[1L] ){
					p <- str_platform(object)
					if( !is.null(p) ) p <- str_c(' - ', p)
				}
				cat("Expression: ", if( res[1L] ) str_c(d[1L], ' x ', d[2L], p) else FALSE, "\n",
						"Components: ", if( d[3L] ) str_c(d[3L], ' [', str_c(selectSome(basisnames(object), 3), collapse=', '), ']') else 0L, "\n",
						"Signatures: ", res[2L], " - Proportions: ", res[3L], "\n", sep='')
			}
			
			invisible(res)
		}
)

#' Subsetting ExpressionMix Objects
#' 
#' Subset method for ExpressionMix objects, which subsets both the 
#' expression and mixture data.
#' 
#' @rdname ExpressionMix-subset
#' @name ExpressionMix-subset
NULL

#' @inheritParams base::[
#' @rdname ExpressionMix-subset
#' @export
setMethod('[', 'ExpressionMix', 
		function (x, i, j, ..., drop = FALSE)
		{	
			if( !missing(i) && missing(j) ){
				x <- selectMethod('[', 'ExpressionSet')(x, i, , drop=drop)
				if( hasBasis(x) )
					x <- selectMethod('[', 'NMF')(x, i, , ..., drop=drop)
			}
			else if( missing(i) && !missing(j) ){
				x <- selectMethod('[', 'ExpressionSet')(x, , j, drop=drop)
				if( hasCoef(x) )
					x <- selectMethod('[', 'NMF')(x, , j, ..., drop=drop)
				slot(x, 'pure') <- slot(x, 'pure')[j] 
			}else if( !missing(i) && !missing(j) ){
				x <- selectMethod('[', 'ExpressionSet')(x, i, j, drop=drop)
				if( hasBasis(x) )
					x <- selectMethod('[', 'NMF')(x, i, , ..., drop=drop)
				if( hasCoef(x) )
					x <- selectMethod('[', 'NMF')(x, , j, ..., drop=drop)
				slot(x, 'pure') <- slot(x, 'pure')[j]
			}else{ # all missing: possibly not missing k for subsetting on basis components			
				x <- selectMethod('[', 'NMF')(x, , , ..., drop=drop)
			}
			
			x
		}
)

#' Creates an ExpressionMix object from an ExpressionSet object.
#' 
#' @param composition an \code{\linkS4class{NMF}} object, that contains the 
#' mixture data (signatures and/or proportions).
#' 
setMethod('ExpressionMix', 'ExpressionSet', 
		function(object, ..., composition=nmfModel(object)){
			
			# create 'empty' ExpressionMix object 
			if( hasBasis(composition) )
				rownames(composition) <- featureNames(object)
			if( hasCoef(composition) )
				colnames(composition) <- sampleNames(object)
			res <- new('ExpressionMix', exprs=exprs(object), ...
					, W=basis(composition), H=coef(composition))
			
			# On has to copy each slot copy slot of ExpressionSet 
			# this is due to the wrong design of initialize,ExpressionSet
			for(s in slotNames(object)){
				slot(res, s) <- slot(object, s)
			}
			
			# pass on any markers attached to the NMF model
			res <- attachMarkers(res, getMarkers(composition))
			
			res
		}
)

#' Creates an ExpressionMix object using the matrix \code{object} as the 
#' expression matrix.
#' @inheritParams ExpressionMix,ExpressionSet-method
setMethod('ExpressionMix', 'matrix', 
		function(object, ..., composition=nmfModel(object)){
			res <- new('ExpressionMix', exprs=object, ...
					, W=basis(composition), H=coef(composition))
			
			# pass on any markers attached to the NMF model
			res <- attachMarkers(res, getMarkers(composition))
			
			res
		}
)

#' Creates an ExpressionMix object from the accession key of a registered dataset.
#' 
#' The expression data is first loaded with \code{\link{gedData}} -- as an
#' ExpressionSet object. The mixture data is then built as an 
#' \code{\linkS4class{NMF}} object and put together with the expression data 
#' into a new \code{\linkS4class{ExpressionMix}} object.  
#' 
#' @param no.pure a logical that indicates if the pure samples should 
#' be filtered out.
#' @param verbose logical that indicates if log messages should be shown.
#' 
setMethod('ExpressionMix', 'character', 
		function(object, no.pure = FALSE, verbose=FALSE, ...){
			
			# set verbosity level
			if( !missing(verbose) ){
				ol <- lverbose(verbose)
				on.exit( lverbose(ol) )
			}
			
			vmessage("Loading dataset '", object, "' ... ", appendLF=FALSE)
			# get the data entry
			d <- gedData(object)		
			# build an 'ExpressionMix' based on this data
			eset <- eset(d, ..., verbose=verbose-1)
			# create pure vector if necessary
			p <- logical()
			if( !is_NA(pidx <- pureIndexes(d)) ){
				p <- rep(FALSE, ncol(eset))
				p[pidx(eset)] <- TRUE
			}
			
			# load signatures
			lmessage(2, "\n Loading basis signatures ... ", appendLF=FALSE)
			b <- basis(d)
			lmessage(2, if( nrow(b) ) 'OK' else 'NONE' 
						, ' [', nrow(b), ' x ', ncol(b), ']')
			if( nrow(b) ) 
				lmessage(3, '  * Signatures: ', str_out(colnames(b)))
			#
			
			# load proportions
			lmessage(2, " Loading proportions ... ", appendLF=FALSE)
			co <- coef(d)
			lmessage(2, if( ncol(co) ) 'OK' else 'NONE'
						, ' [', nrow(co), ' x ', ncol(co), ']')
			if( ncol(co) )
				lmessage(3, '  * Proportions: ', str_out(rownames(co)))
			#
			
			# make object
			lmessage(2, " Warpping into an ExpressionMix object ... ", appendLF=FALSE)
			res <- ExpressionMix(eset, pure=p, composition=nmfModel(b, co))
			lmessage(2, "OK")
			
			# filter out pure samples
			if( no.pure ) res <- mixedSamples(res)
			vmessage("OK")
			
			res
		}
)

#' Extracts an ExpressionSet object from an \code{ExpressionMix} object.
#' 
#' Because \code{\linkS4class{ExpressionMix}} inherits from 
#' \code{\linkS4class{ExpressionSet}}, this methods simply returns the 
#' object coerced as an \code{ExpressionSet} object.
#' 
#' @examples
#' 
#' \dontrun{
#' # get data from the registry
#' e <- gedData('GSE20300')
#' # this is an ExpressionMix object
#' e
#' # extract the ExpressionSet part
#' class(eset(e))
#' }
#' 
setMethod('eset', 'ExpressionMix', 
	function(object, ...){
		as(object, 'ExpressionSet')
	}
)

#' Extracts mixture data from an \code{ExpressionMix} object.
#' 
#' Because \code{\linkS4class{ExpressionMix}} inherits from 
#' \code{linkS4class{NMFstd}}, this methods simply returns the
#' the object coerced as an \code{NMFstd} object.
#' 
#' @examples
#' 
#' e <- ExpressionMix(rmatrix(10,5))
#' e
#' # the mixture data is stored as an NMF model
#' mixData(e)
#' 
#' # assign a random mixture data with 5 cell types
#' m <- rnmf(3, e)
#' # ensure proportions sum-up to one
#' coef(m) <- scoef(m)
#' mixData(e) <- m
#' 
setMethod('mixData', 'ExpressionMix', 
	function(object, ...){
		as(object, 'NMFstd')
	}
)

#' Sets mixture data in an ExpressionMix object.
#' 
#' Because \code{\linkS4class{ExpressionMix}} inherits from 
#' \code{linkS4class{NMFstd}}, this method simply copies all 
#' slots in \code{value} in into the corresponding slots in 
#' \code{object}, and then checks that the result object is valid. 
#' 
setReplaceMethod('mixData', signature('ExpressionMix','NMFstd'), 
	function(object, value){
		# copy all slots from value
		for(s in slotNames(value)){
			slot(object, s) <- slot(value, s)
		}
		validObject(object)
		object
	}
)

#' Subsetting Data with MarkerList Objects
#'
#' The function \code{subsetML} enables direct subsetting of data objects with 
#' \code{\link{MarkerList}} objects, as long as they have a suitable
#' \code{\link{featureNames}}, \code{\link{rownames}} or \code{\link{names}} 
#' method, as well as a subset method \code{[}.
#' It is used for defining \code{'['} methods for \code{\linkS4class{NMF}} 
#' and \code{\linkS4class{ExpressionMix}}.
#'  
#' @inheritParams base::[
#' @param drop logical passed to the subset method \code{'['} of \code{x}.  
#' 
#' @export
#' @examples
#' 
#' data(sample.ExpressionSet)
#' e <- sample.ExpressionSet
#' annotation(e)
#' 
#' 
subsetML <- function (x, i, j, ..., drop = FALSE){
	
	# unlist ids
	ids <- marknames(i)
	
	if( is.integer(ids) ){ # direct subset with index
		ogenes <- ids
		
	}else{ # find name method to use if necessary
		
		# convert if necessary
		i <- tryConvertIDs(i, x)
		ids <- marknames(i)
		#
		
		cl <- class(x)
		if( is.null(nfun <- selectMethod('featureNames', cl, optional=TRUE)) ){
			if( is.null(nfun <- selectMethod('rownames', cl, optional=TRUE)) ){
				if( is.null(nfun <- selectMethod('names', cl, optional=TRUE)) ){
					stop("Could not find a suitable method `featureNames`, `rownames` or `names` for object of class '", cl, "'.")
				}
			}
		}
		ogenes <- ids[ids %in% nfun(x)]
		
	}		

  # it is important not to pass drop if missing as some subset methods rely on the number of arguments
	if (missing(j)){
		if( missing(drop) ) x[ogenes, ,...]
		else x[ogenes, , ..., drop=drop]
	}else{
		if( missing(drop) ) x[ogenes, j, ...]
		else x[ogenes, j, ..., drop=drop]
	}
}

#' @rdname ExpressionMix-subset
#' @export
setMethod('[', signature(x="MatrixData", i="MarkerList", j="ANY", drop="ANY")
	, subsetML
)
#' @rdname ExpressionMix-subset
#' @export
setMethod('[', signature(x="NMF", i="MarkerList", j="ANY", drop="ANY")
	, subsetML
)

#' Sets the feature names on both the \code{ExpressionSet} and 
#' \code{NMF} objects.
#' @export
#' @rdname ExpressionMix-dims
setReplaceMethod('featureNames', "ExpressionMix", 
	function(object, value){
		rownames(object) <- value
		object
	}
)

#' Sets the sample names on both the \code{ExpressionSet} and 
#' \code{NMF} objects.
#' @export
#' @rdname ExpressionMix-dims
#' @aliases sampleNames<-,ExpressionMix,ANY-method
setReplaceMethod('sampleNames', "ExpressionMix", 
	function(object, value){
		colnames(object) <- value
		object
	}
)
