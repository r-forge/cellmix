# TODO: Add comment
# 
# Author: Renaud Gaujoux
# Created: Mar 6, 2013
###############################################################################

#' @include registry-algorithms.R
NULL

#' Partial Gene Expression Deconvolution: Marker Mean Expression Profile
#' 
#' @description
#' The algorithm \sQuote{meanProfile} uses a set of known marker genes 
#' for each cell type to compute a mean expression profile within each cell type 
#' separately.
#' 
#' This is in essence the preliminary step proposed by \cite{Kuhn2011} in order to 
#' compute a proxy for the actual cell proportions.
#' 
#' The average profiles are expected to be correlate well with the actual proportions,
#' provided that the individual markers gene expression profiles are not 
#' too noisy and that the markers are indeed markers.
#' 
#' @details
#' \strong{Important:} this method does not compute cell-specific differential expression
#' as described in \cite{Kuhn2011}, but only the cell proportion proxy.
#' Hence, the result only contains estimated proportions (accessible via \code{\link{coef}}) 
#' and an empty basis signature matrix. 
#' 
#' @param y target expression matrix
#' @param x specification of which feature to use as markers, as a list,
#' typically a \code{\link{MarkerList}} object.
#' @param scale logical that indicates if the proportion estimates should 
#' be scale to sum-up to one before returning them.
#' @param ... extra arguments not used.
#' @aliases meanProfile-ged
#' @examples
#' 
#' # random data with markers 
#' x <- rmix(3)
#' m <- getMarkers(x)
#' 
#' # compute proxy proportions
#' res <- ged(x, m, method='meanProfile')
#' # no cell-specific signatures
#' dim(res)
#' 
#' #NB: estimates are not scaled to sum up to one
#' profplot(x, res)
#' profplot(x, res, scale=TRUE)
#' 
gedAlgorithm.meanProfile <- setGEDMethod(key='meanProfile'
		, description = "Compute proportion proxies as mean expression profiles"
		, algorithm = function(y, x, scale=FALSE, ...){
			
			data <- x
			i <- matchIndex(data, y)
			# compute proportions
			if( nmark(i) ){
				p <- colMeansBy(exprs(y), i)
				if( scale ) p <- scoef(p)
				p
			}else{
				warning("Could not match any of the markers [", str_out(marknames(data)), ']'
						, " to features in the target data [", featureNames(y), ']')
				NA_Matrix(length(data), 0)
			}
		}
		, reqBasis = FALSE
		, reqCoef= FALSE
		, reqMarker = TRUE
		, maxIter = 1L
)

