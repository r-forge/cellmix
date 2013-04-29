# Stats utility functions
# 
# Author: Renaud Gaujoux
# Creation: 23 Jan 2012
###############################################################################

#library(matrixStats)

.applyBy_BY <- pkgmaker::sVariable(NULL)

#' Group Apply
#' 
#' \code{appplyBy} is an S3 generic function that applies a given function 
#' to sub-matrices of a matrix-like object, which are generated according to 
#' a factor that defines groups rows or columns.
#' 
#' @export
applyBy <- function(x, ...){
	UseMethod('applyBy')
}

#' The method \code{applyBy.matrix} is the work horse function 
#' that is called by other more user-friendly functions.
#' 
#' \code{applyBy.matrix} is a wrapper around \code{\link[matrixStats]{colAvgsPerRowSet}}, 
#' which make the computation really fast, but requires somehow cumbersome matrix 
#' specifications for the groups of columns or rows.
#' The wrapper builds the arguments for the particular case where 
#' the groups are defined by a factor.
#' 
#' @param x matrix-like object on which \code{\link{apply}} can be called.
#' @param BY factor or object coerced to a factor, that defines the groups within 
#' which the function \code{FUN} is applied.
#' 
#' If \code{x} is an ExpressionSet object, then \code{BY} can be the names of a
#' sample (resp. feature) annotation variable if \code{MARGIN=1} (resp. \code{MARGIN=2L}) 
#' (see examples).
#' @param MARGIN margin along which the function \code{FUN} is applied: 
#' 1L for rows, 2L for columns. 
#' @param FUN function to apply to each sub-matrix that contains the rows/columns 
#' defined by each level of argument \code{BY}.
#' It must be a function that takes a matrix as its first argument and returns a vector 
#' of length the dimension of margin \code{MARGIN} of \code{x}.
#' @inheritParams matrixStats::colAvgsPerRowSet
#' @param ... extra parameters passed to \code{FUN}.
#' @param DROP logical that indicates if absent levels should be removed 
#' from the result matrix, or appear as 0-filled rows/columns.
#' 
#' @return The result is a matrix or an \code{ExpressionSet} object 
#' whose margin's dimension \code{MARGIN} is equal the same margin's 
#' dimension in \code{x}, and the other to the number of levels 
#' in \code{BY}.
#'
#' @S3method applyBy matrix 
#' @importFrom matrixStats colAvgsPerRowSet rowAvgsPerColSet
#' @importFrom matrixStats colAvgsPerRowSet.matrix rowAvgsPerColSet.matrix
#' @rdname applyBy 
#' @examples
#' 
#' # random data matrix
#' x <- rmatrix(12, 6)
#' 
#' # by groups of columns
#' fc <- gl(2, 3)
#' b <- applyBy(x, fc, 1L, rowSums)
#' b
#' # or
#' balt <- rowApplyBy(x, fc, rowSums)
#' identical(b, balt)
#' 
#' # by groups of rows
#' fr <- gl(3, 4)
#' b <- applyBy(x, fr, 2L, colSums)
#' # or
#' balt <- colApplyBy(x, fr, colSums)
#' identical(b, balt)
#'  
applyBy.matrix <- function(x, BY, MARGIN, FUN, W=NULL, ..., DROP=FALSE){
	
	BYMARGIN <- 3L-MARGIN
	# check arguments
	mdim <- length(dim(x))
	if( MARGIN > mdim )
		stop("Invalid value for argument `MARGIN`: greater than the number of dimensions of `x` [", mdim, "]")	
	if( BYMARGIN > mdim )
		stop("Invalid value for argument `BYMARGIN`: greater than the number of dimensions of `x` [", mdim, ']')
	if( MARGIN > 2 || BYMARGIN > 2 )
		stop("Invalid margins: must be either 1L or 2L [", MARGIN, ', ', BYMARGIN, ']')
	
	# convert list of indexes into a vector/factor
	if( is.list(BY) ){ 
		idx <- unlist(BY)
		if( !is.integer(idx) || (is.character(idx) && any(!idx %in% rownames(x))) )
			stop("Invalid `BY` argument of type list [", str_out(idx), "]:"
				, " should contain integer indexes or character strings that match rownames in the data [", str_out(rownames(x)), "].")
		if( anyDuplicated(idx) )
			stop("Invalid `BY` argument of typ list: it should not contain any duplicated values.")
		# use names as levels
		if( is.null(names(BY)) ) names(BY) <- seq_along(BY)		
		# subset to the given indexes
		x <- if( BYMARGIN == 1L ) x[idx, ,drop=FALSE] else x[, idx, drop=FALSE]
		# convert into a factor
		v <- unlist(mapply(rep, names(BY), sapply(BY, length)))
		BY <- factor(v, levels=names(BY))
	}
	
	bydim <- dim(x)[BYMARGIN]
	if( length(BY) != bydim )
		stop("Invalid value for argument `BY`: length [",length(BY),"] is not equal to the dimension [", bydim, "] of the BY margin [", BYMARGIN, ']')
	
	# coerce to factor if necessary
	if( !is.factor(BY) ) BY <- factor(BY, levels=unique(BY))	
	
	# build subset matrix
	s <- split(1:bydim, BY)
	nm <- max(sapply(s, length))
	S <- matrix(0, nm, length(s))
	colnames(S) <- names(s)
	sapply(seq_along(s), function(i){
		idx <- s[[i]]
		if( length(idx) > 0L ) S[1:length(idx),i] <<- idx
	})
	# store idx in static variable if requested
	if( isTRUE(.applyBy_BY()) ){
		.applyBy_BY(BY)
	}
	
	# call relevant function from matrixStats
  res <- 
	if( MARGIN == 1L ) rowAvgsPerColSet(X=x, S=S, FUN=FUN, W=W, ...)
	else colAvgsPerRowSet(X=x, S=S, FUN=FUN, W=W, ...)

  # drop absent levels if requested
  if( DROP ){
	  lv <- levels(BY)
	  if( length(w <- which(!lv %in% unique(as.character(BY)))) ){
		# update stored factor if necessary
		if( !is.null(.applyBy_BY()) ){
		  .applyBy_BY(droplevels(BY))
		}
		if( MARGIN == 1L ){
			res <- res[,-w,drop=FALSE]
#			tmp <- cbind(res, NA_Matrix(nrow(res), length(w)))
#			colnames(tmp)[(ncol(res)+1):ncol(tmp)] <- lv[w] 
		}else{
			res <- res[-w,,drop=FALSE]
#			tmp <- rbind(res, NA_Matrix(length(w), ncol(res)))
#			rownames(tmp)[(nrow(res)+1):nrow(tmp)] <- lv[w]
		   }
		}
	}

	#
	
	# return result
	res
}

#' A method is provided for \code{\link{ExpressionSet}} objects, 
#' which preserve sample and feature annotations.
#' Moreover it allows directly passing names of feature/sample annotation -- factor -- variables 
#' in argument \code{BY} (see examples).
#' 
#' @param ANNOTATIONS logical that indicates if samples/feature annotations should 
#' be kept, when the input data is an \code{\link{ExpressionSet}} object.
#' Currently, if \code{TRUE}:
#' \itemize{
#' \item if code{MARGIN=1L}, then feature annotations are kept unchanged, and 
#' phenotypic sample annotations are discarded. 
#' \item if code{MARGIN=2L}, then phenotypic sample annotations are kept unchanged, and 
#' feature annotations are discarded.
#' } 
#' 
#' In any case, the value of slot \code{annotation} (i.e. the annotation package), 
#' is passed on to the result object.
#'  
#' @S3method applyBy ExpressionSet 
#' @rdname applyBy
#' @examples
#' 
#' ## Method for ExpressionSet objects
#' 
#' x <- ExpressionSet(x, annotation='abcd.db')
#' y <- rowMinsBy(x, fc)
#' \dontshow{ 
#' 	stopifnot( is(y, 'ExpressionSet') )
#'  stopifnot( identical(nrow(y), nrow(x)) )
#'  stopifnot( ncol(y) == nlevels(fc) )
#'  stopifnot( identical(annotation(y), 'abcd.db') )
#' }
#' y <- colMinsBy(x, fr)
#' \dontshow{ 
#' 	stopifnot( is(y, 'ExpressionSet') )
#'  stopifnot( identical(ncol(y), ncol(x)) )
#'  stopifnot( nrow(y) == nlevels(fr) ) 
#'  stopifnot( identical(annotation(y), 'abcd.db') )
#' }
#' 
#' ## annotations are conserved/collapsed
#' pData(x) <- data.frame(Group=fc, Sample=letters[1:ncol(x)])
#' pData(x)
#' fData(x) <- data.frame(ENTREZID=fr, Gene=letters[nrow(x):1])
#' fData(x)
#' 
#' # keep feature annotations, collapse sample annotations
#' y <- rowMinsBy(x, 'Group')
#' pData(y)
#' fData(y)
#' \dontshow{ 
#' 	stopifnot( is(y, 'ExpressionSet') ) 
#'  stopifnot( identical(nrow(y), nrow(x)) )
#'  stopifnot( ncol(y) == nlevels(fc) )
#'  stopifnot( identical(annotation(y), 'abcd.db') )
#'  stopifnot( identical(fData(y), fData(x)) )
#'  stopifnot( nrow(pData(y)) == nlevels(fc) )
#' }
#' 
#' # keep sample annotations, collapse feature annotations 
#' y <- colMinsBy(x, 'ENTREZID')
#' pData(y)
#' fData(y)
#' \dontshow{ 
#' 	stopifnot( is(y, 'ExpressionSet') ) 
#'  stopifnot( identical(annotation(y), 'abcd.db') )
#'  stopifnot( identical(ncol(y), ncol(x)) )
#'  stopifnot( nrow(y) == nlevels(fr) ) 
#'  stopifnot( identical(pData(y), pData(x)) )
#'  stopifnot( nrow(fData(y)) == nlevels(fr) )
#' }
applyBy.ExpressionSet <- function(x, BY, MARGIN, ..., ANNOTATIONS=TRUE){
	
	# convert single character string into annotation variable
	if( isString(BY) ){
		if( MARGIN == 1L ){ # phenotypic variable
			if( !BY %in% varLabels(x) ){
				stop("Invalid string argument BY: there is no phenotypic/sample annotation variable called '", BY, "'.\n"
					, "  Defined variable are: ", str_out(varLabels(x), Inf))
			}
			BY <- pData(x)[[BY]]
		}else{ # feature annotation variable
			if( !BY %in% names(fData(x)) ){
				stop("Invalid string argument BY: there is no feature annotation variable called '", BY, "'.\n"
						, "  Defined variable are: ", str_out(names(fData(x)), Inf))
			}
			BY <- fData(x)[[BY]]
		}
	}
		
	# apply to expression matrix
	.applyBy_BY(TRUE)
	on.exit( .applyBy_BY(NULL) )
	res <- applyBy(exprs(x), BY=BY, MARGIN=MARGIN, ...)
	
	# re-wrap into an ExpressionSet objects
	library(Biobase)
	# pass on annotations whenever possible
	fd <- pd <- NULL
	if( ANNOTATIONS ){
		if( MARGIN == 1L ){
			if( nrow(ad <- featureData(x)) > 0L ) fd <- ad # keep feature annotations
			# collapse sample annotation for non-list BY arguments
			if( !is.list(BY) && nrow(ad <- phenoData(x)) > 0L ){
				# get used BY factor
				fBY <- .applyBy_BY()
				# extract annotation for first representative of each level
				df <- pData(x)[!duplicated(as.character(fBY)) & !is.na(fBY), ]
				rownames(df) <- sampleNames(res)
				pd <- AnnotatedDataFrame(df)
			}
		}else if( MARGIN == 2L ){
			if( nrow(ad <- phenoData(x)) > 0L ) pd <- ad # keep sample annotations
			# collapse feature annotation for non-list BY arguments
			if( !is.list(BY) && nrow(ad <- featureData(x)) > 0L ){
				# get used BY factor
				fBY <- .applyBy_BY()
				# extract annotation for first representative of each level
				df <- fData(x)[!duplicated(as.character(fBY)) & !is.na(fBY), ]
				rownames(df) <- featureNames(res)
				fd <- AnnotatedDataFrame(df)
			}
		}
	}
	# no sample/feature annotation
	ca <- call('ExpressionSet', res, annotation=annotation(x))
	if( !is.null(pd) ) ca$phenoData <- pd
	if( !is.null(fd) ) ca$featureData <- fd
	res <- eval(ca)
	
	res
}

#' \code{rowApplyBy} applies a function to rows of sub-matrices whose columns 
#' are defined by a factor.
#' @export
#' @rdname applyBy
rowApplyBy <- function(x, BY, FUN, ...){ 
	applyBy(x, BY=BY, MARGIN=1L, FUN=FUN, ...)
}

#' \code{rowApplyBy} applies a function to columns of sub-matrices whose rows 
#' are defined by a factor.
#' @export
#' @rdname applyBy
colApplyBy <- function(x, BY, FUN, ...){ 
	applyBy(x, BY=BY, MARGIN=2L, FUN=FUN, ...)
}



# generator functions
.rowApplyByFunction <- function(FUN){
  function(x, BY, ...){
    applyBy(x, BY=BY, MARGIN=1L, FUN=FUN, ...)
  }
}
.colApplyByFunction <- function(FUN){
  function(x, BY, ...){
    applyBy(x, BY=BY, MARGIN=2L, FUN=FUN, ...)
  }
}

#' \code{col<STAT>By} computes for each column a given statistic within separate groups of rows, which are defined by a factor.
#' @export
#' @rdname applyBy
colSumsBy <- .colApplyByFunction(colSums)

#' \code{row<STAT>By} computes for each row a given statistic within separate groups of columns, which are defined by a factor.
#' @export
#' @rdname applyBy
rowSumsBy <- .rowApplyByFunction(rowSums)

#' @export
#' @rdname applyBy
rowMeansBy <- .rowApplyByFunction(rowMeans)
#' @export
#' @rdname applyBy
colMeansBy <- .colApplyByFunction(colMeans)

#' @export
#' @rdname applyBy
rowMediansBy <- .rowApplyByFunction(rowMedians)
#' @export
#' @rdname applyBy
colMediansBy <- .colApplyByFunction(colMedians)

#' @export
#' @rdname applyBy
rowMaxsBy <- .rowApplyByFunction(rowMaxs)
#' @export
#' @rdname applyBy
colMaxsBy <- .colApplyByFunction(colMaxs)

#' @export
#' @rdname applyBy
rowMinsBy <- .rowApplyByFunction(matrixStats::rowMins)
#' @export
#' @rdname applyBy
colMinsBy <- .colApplyByFunction(matrixStats::colMins)
