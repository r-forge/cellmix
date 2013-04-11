# Stats utility functions
# 
# Author: Renaud Gaujoux
# Creation: 23 Jan 2012
###############################################################################

library(matrixStats)

#' Group Apply
#' 
#' Applies a function to rows or columns of a matrix-like object, separately 
#' within groups defined by a factor.
#' 
#' \code{applyBy} is the work horse function that is called by other more 
#' user-friendly functions.
#' 
#' @details
#' \code{applyBy} is a wrapper around \code{\link[matrixStats]{colAvgsPerRowSet}}, 
#' which make the computation really fast, but requires somehow cumbersome matrix 
#' specifications for the groups of columns or rows.
#' The wrapper builds the arguments for the particular case where 
#' the groups are defined by a factor.
#' 
#' @param x matrix-like object on which \code{\link{apply}} can be called.
#' @param BY factor or object coerced to a factor, that defines the groups within 
#' which the function \code{FUN} is applied.
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
#' @return The result is a matrix whose margin's dimension \code{MARGIN} is equal 
#' the same margin's dimension in \code{x}, and the other to the number of levels 
#' in \code{BY}. 
#' 
#' @export
#' 
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
applyBy <- function(x, BY, MARGIN, FUN, W=NULL, ..., DROP=FALSE){
	
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
		x <- x[idx, ]
		# convert into a factor
		v <- unlist(mapply(rep, names(BY), sapply(BY, length)))
		BY <- factor(v, levels=names(BY))
	}
	
	bydim <- dim(x)[BYMARGIN]
	if( length(BY) != bydim )
		stop("Invalid value for argument `by`: length [",length(BY),"] is not equal to the dimension [", bydim, "] of the by margin [", BYMARGIN, ']')
	
	# coerce to factor if necessary
	if( !is.factor(BY) ) BY <- factor(BY, levels=unique(BY))	
	
	# special case for ExpressionSet object
  	x.orig <- x
	if( is(x, 'ExpressionSet') ){
    	x <- exprs(x)
	}
  
	# build subset matrix
	s <- split(1:dim(x)[3-MARGIN], BY)
	nm <- max(sapply(s, length))
	S <- matrix(0, nm, length(s))
	colnames(S) <- names(s)
	sapply(seq_along(s), function(i){
		idx <- s[[i]]
		if( length(idx) > 0L ) S[1:length(idx),i] <<- idx
	})
	
	# call relevant function from matrixStats
  res <- 
	if( MARGIN == 1L ) rowAvgsPerColSet(X=x, S=S, FUN=FUN, W=W, ...)
	else colAvgsPerRowSet(X=x, S=S, FUN=FUN, W=W, ...)

  # complete absent levels with NAs if requested
  if( DROP ){
	  lv <- levels(BY)
	  if( length(w <- which(!lv %in% unique(as.character(BY)))) ){
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
  
  # re-wrap ExpressionSet objects
	if( is(x.orig, 'ExpressionSet') ){
	  res <- ExpressionSet(res, annotation=annotation(x.orig))
	}
  
  # return result
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
