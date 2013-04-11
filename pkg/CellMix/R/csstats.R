# Cell-specific statistics
# 
# Author: Renaud Gaujoux
# Created: 11 Dec 2012
###############################################################################


#' Compute Cell-Specific Statistics
#' 
#' \code{csTopTable} is a generic function that returns statistics 
#' for each genes, at the cell type levels.
#' It is used on deconvolution results from \code{\link{ged}}. 
#' 
#' @param x data object, typically returned by \code{\link{ged}}.
#' @param ... extra parameters passed to subsequent calls.
#' 
#' @export
csTopTable <- function(x, ...){
	UseMethod('csTopTable')
}

dimidx <- function(x, margin){
	d <- dimnames(x)[[margin]]
	if( is.null(d) ) d <- seq(1, dim(x)[[margin]]) 
	d
}

#' @param n maximum number of features to extract -- in each cell type.
#' @param ct specifies the cell type for which one wants to 
#' extract the top features.
#' @param decreasing logical that indicates the feature ordering, based on 
#' their p-values or FDRs. 
#' 
#' @S3method csTopTable matrix
#' @rdname csTopTable
csTopTable.matrix <- function(x, n=100L, ct=NULL, decreasing=FALSE, ...){
	
	cts <- dimidx(x, 2L)
	top <- sapply(cts, function(i){
		p <- setNames(x[, i], dimidx(x, 1L))
		p <- p[order(p, decreasing=decreasing)]
		if( !is.null(n) ) p <- head(p, n)
		p
	}, simplify=FALSE)

	# select single cell type
	if( !is.null(ct) ) top <- top[ct]
	if( length(top) == 1L ) top <- top[[1L]] 
	
	top
}

#' @S3method csTopTable array
#' @rdname csTopTable
csTopTable.array <- function(x, ...){
	
	if( length(dim(x)) == 2L ) csTopTable(x[,], ...)
	else{
		cmp <- dimidx(x, 1L)
		sapply(cmp, function(i, ...){
			csTopTable(x[i,,], ...)
		}, ..., simplify=FALSE)
	}
}

#' @S3method csTopTable NMFfit
#' @rdname csTopTable
csTopTable.NMFfit <- function(x, ...){
	x <- structure(list(x), class=str_c('ged_', algorithm(x)))
	csTopTable(x, ...)
}

#' @S3method csTopTable NMFfitX
#' @rdname csTopTable
csTopTable.NMFfitX <- function(x, ...){
	csTopTable(minfit(x), ...)
}
