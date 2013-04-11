# Assigning components to cell types
# 
# Author: Renaud Gaujoux
# Created: 29 Nov 2012
###############################################################################


#' Assigning Estimated Signature to Real Cell-Types
#' 
#' Assigns and reorder basis components to a cell type based on the 
#' proportion of cell-specific markers that contribute the most to each component.
#' 
#' @param object an \code{NMF} object
#' @param reference a list of markers, usually a \code{MarkerList} object.
#' @param n the number of markers to use per cell-type to actually assign the components
#' @return \code{object} with the basis components reordered -- and renamed -- according to the 
#' computed assigned cell-types.
#' 
#' @export
match.nmf <- function(object, reference, n=NA){
	
	# assign each component of object to one in the reference
	t <- match.type(object, reference=reference, n=n, return.table=TRUE)
	idx <- t$match
	# rename the basis
	basisnames(object) <- colnames(t$table)[idx]
	# reverse the mapping
	idx <- match(seq_along(idx), idx, 0L)
	if( length(idx) != nbasis(object) )
		warning("Not all the components were assigned to a type.")
	# reorder the original object
	object[,,idx]
}

#vector that contains the indexes of the true cell types, one for each estimated component, so that
#if \code{res <- match.type(x, ref)} then x[,res
setGeneric('match.type', function(object, ...) standardGeneric('match.type'))
setMethod('match.type', 'NMF'
	, function(object, ...){
		match.type(basis(object), ...)		
	}
)
setMethod('match.type', 'list'
	, function(object, reference){
		res <- match.type(reference, object, ..., return.table=TRUE)
		match(1:ncol(res$table), res$match)
	}
)
setMethod('match.type', 'matrix'
	, function(object, reference, n, verbose=FALSE, return.table=FALSE, use.percentage=TRUE, ...){
		#, function(object, reference, n, threshold=0, verbose=FALSE, return.table=FALSE, use.percentage=TRUE, ...){
		
		# extract the list of markers from the reference object 
		marks <- MarkerList(reference)			
		
		if( length(marks) != ncol(object) )
			stop("match.type - Incompatible dimensions: the number of marker type (", length(marks), ") should be equal to the number of basis (", ncol(object), ")")
		
		# subset the markers if required
		if( !missing(n) && !is.na(n) )
			marks <- subset(marks, n)
		
		m.idx <- matchIndex(marks, object)	
		# compute the true and estimated mapping 
		mark.map <- unlist(lapply(seq_along(marks), function(i) rep(i, length(m.idx[[i]]))))
		m.idx <- unlist(m.idx)
		# the mappings should have the same length
		if( length(m.idx) != nmark(marks) )
			warning("match.type - Some of the reference markers are not present in the data.")
		mark.map.estim <- max.col(object[m.idx,])
		#cat('True :', mark.map, "\n"); cat('Estim:', mark.map.estim, "\n");
		
		# compute the contingency table: estimate in rows, reference in columns
		tmark <- table( estim=factor(mark.map.estim, levels=1:ncol(object))
				, ref=factor(mark.map, levels=seq_along(marks)))		
		
		# use percentages instead of raw counts
		if( use.percentage ){
			n.mark <- sapply(marks, length)
			tmark <- sweep(tmark, 2, 100 / ifelse(n.mark==0, 1, n.mark), '*')
		}
		
		# add names if present
		if( !is.null(names(marks)) )
			colnames(tmark) <- names(marks)
		if( !is.null(colnames(object)) )
			rownames(tmark) <- colnames(object)
		
		if( verbose ){
			cat("Marker distribution:\n")
			print(round(tmark,2))
			#cat("Using threshold =", threshold, "\n\n")
		}
		
		# for each estimated component look for the major class of markers		
		T.tmp <- tmark				
		res <- rep(0, ncol(T.tmp))
		for( i in 1:ncol(T.tmp) ){
			# get the row and column index of the maximum over the remaining entries 
			xm <- which.max(T.tmp)-1
			jm <- xm %/% nrow(T.tmp) + 1
			im <- xm - (jm-1) * nrow(T.tmp) + 1
			
			# assign the estimate row to the inferred reference column
			stopifnot( res[im]==0 )
			res[im] <- jm
			
			# erase the assigned estimate row
			T.tmp[im,] <- NA
			# erase the assigned reference column
			T.tmp[,jm] <- NA
		}
		
		# if there is more basis than cell-type then assign the remaining components
		# to their most contributing type.
#		if( length(reference) != ncol(object) )
#			NA;
		
		# return the mapping as an integer vector
		res <- as.integer(res)
		if( return.table )
			res <- list(match=res, table=tmark)
		
		# return result
		res
	}
)

#' \code{match.cell} is an alias for \code{match.nmf}
#' @export
#' @rdname match.nmf 
match.cell <- match.nmf
