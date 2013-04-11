# Plotting functions
# 
# Author: Renaud Gaujoux
# Creation: 12 Jan 2012
###############################################################################

#' @include AllClasses.R
#' @include markers.R
#' @include ExpressionMix-class.R
NULL

#' Plotting Markers
#' 
#' The method \code{beeswarm} for \code{MarkerList} objects draws, 
#' separately for each cell-type, a stripchart of the markers' 
#' expression values in a reference data, highlighting the expression 
#' of their respective associated -- generally pure -- samples.
#'  
#' @param x an object of class \code{\link{MarkerList}}
#' 
#' @param data a numeric \code{matrix} or an \code{ExpressionSet} object 
#' containing expression values from samples.
#'  
#' @param types a \code{factor} giving the cell-type of each sample
#' @param tcol type colour, used for the marker values in samples of 
#' their corresponding types.
#' @param otcol out-of-type colour, used for maker values in samples 
#' that are not from their corresponding type.
#' @param las specifies the orientation of the y-labels: 1 for horizontal, 
#' 2 for vertical. See \code{\link{par}}.
#' @param fold logical that indicates if the plot should show the markers' plain 
#' expression values or their fold change against the maximum expression 
#' value in samples from other types.
#' 
#' @inheritParams beeswarm::beeswarm
#' 
#' @param ... extra parameters passed to subsequent calls to \code{\link{beeswarm}}, 
#' \code{\link{matplot}}, \code{\link{barplot}}, \code{\link{hist}} or 
#' \code{\link{plot}}, depending on the main calling function.
#' 
#' @rdname markerPlots
#' 
#' @import beeswarm
#' @S3method beeswarm MarkerList
beeswarm.MarkerList <- function(x, data, types=NULL, tcol='red', otcol='#00000090', las=2, fold=FALSE
							, method='square', corral='wrap', pch=19, pwcol=NULL
							, ylab = if(fold) 'Fold changes' else 'Expression values'
							, ...){
	
	# extract expression values
	data <- exprs(data)
	# match markers on data
	i <- matchIndex(x, data)
	if( !nmark(i) ){
		warning("Data does not contain any of the markers from the list.")
		return(invisible())
	}
	
	# infer type if possible
	if( is.null(types) ){
		if( is.null(cn <- colnames(data)) ){
			if( length(x) != ncol(data) ){
				stop("Could not infer types: list and data have no names, and their lengths are incompatible."
					, "\n  NOTE: number of marker type [", length(x),"] must be equal the number of data columns [", ncol(data),"].")
		}
			types <- names(x)
		}else if( all(cn %in% names(x)) ){
			types <- colnames(types)	
		}else if( any(cn %in% names(x)) ){
			types <- ifelse(cn %in% names(x), cn, NA)
		}else if( !any(cn %in% names(x)) ){
			stop("Could not infer sample types from data column names:"
				, " none of the data column (", str_out(cn), ") match any of the marker types (", str_out(names(x)), ").")
		} 
	}
	
	if( is.null(types) ){
		stop("Could not infer sample types: please provide argument `types`.")
	}
	
	# convert into a factor
	types <- rep(types, length.out=ncol(data))
	if( !is.factor(types) ) types <- factor(types, levels=unique(types))
	# check consistency of levels
	if( !all( levels(types) %in% names(x) ) ){
		types[!types %in% names(x)] <- NA
		if( pb <- sum(is.na(types)) ){
			if( pb == length(types) ){
				stop("Inconsistent sample types: marker list contains no markers for any of the sample types ", str_out(levels(types)), "."
					, "\n  Available types are: ", str_out(names(x)))
			}
		}
	}
	
	# convert factor into list
	ltypes <- split(seq_along(types), types, drop=TRUE)
	nt <- length(ltypes)
	opar <- par(mfrow=c(sqrt(nt),ceiling(sqrt(nt))))
	on.exit(par(opar))
	# subset markers to the ones that have matching samples
	i <- i[names(ltypes), ]
	
	res <- sapply(names(i), 
			function(k){
				intype <- ltypes[[k]]
				idx <- i[[k]]
				n <- length(idx)
				ldata <- lapply(idx, function(i){
					if( !fold ) data[i, ]
					else data[i, ] / max(data[i, -intype])
				})
				if( is.null(pwcol) ){
					pwcol <- rep(otcol, ncol(data))
					pwcol[intype] <- tcol
					pwcol <- replicate(n, pwcol, simplify=FALSE)
				}
				beeswarm(ldata, pwcol = pwcol
						, method = method, corral = corral, pch = pch
						, las=las, main = paste(k, ' (', n, ')', sep='')
						, ylab = ylab
						, ...)
			})
	
	invisible(res)
}

#' The method \code{stripchart.MarkerList} is an alias for \code{beeswarm.MarkerList}. 
#' 
#' @rdname markerPlots
#' @S3method stripchart MarkerList
stripchart.MarkerList <- beeswarm.MarkerList

#' \code{profplot} plots the expression profiles of markers in the list. 
#' 
#' \code{profplot} is vectorised over arguments \var{scale} and \var{legend}.
#'  
#' @inheritParams NMF::profplot
#' @param y reference data matrix, that contains the marker expression values 
#' that are plotted.
#' @param groups factor coercible vector that defines groups of samples, 
#' which are plotted using different point types. 
#' @param scale when \code{split=FALSE}, this is a logical that indicates if 
#' the marker expression profiles should be scaled into relative contribution 
#' (i.e., sum up to one), as in \code{\link{profplot}}.
#' 
#' When \var{split} is not \code{FALSE}, then \code{scale=TRUE} indicates that 
#' the profiles should be normalised with the mean expression profile of all markers.
#' If a single numeric, then it indicates the index of a specific marker 
#' that is used to normalise the expression profiles of other markers.
#' Separate indexes for each cell type can be passed as a numeric vector or a list.
#' 
#' This is useful to highlight groups of markers that should have equivalent 
#' cell-specific expression levels, given the -- mixed -- expression data in \code{y}.
#' These groups are likely to provide better proportion/cell-specific estimates.
#' 
#' Note: values are recycled if necessary. 
#' @param col colour specification for each cell types if \code{split=FALSE}, 
#' or for each marker in each cell type otherwise. 
#' @param split logical that indicates if one should draw a plot for each type 
#' separately.
#' 
#' One can specify how to split the plots, by providing a 2-long vector that is
#' used as \code{par(mfrow=split)}.
#' @param ylab label for the y-axis.
#' If missing, a default label is generated, which, in particular, tells if the plot 
#' is scaled or not.
#' 
#' @param restore.gpar logical that indicates if all graphical parameters should be restored 
#' on exit. 
#' It could be used to add things to the plots, especially when \code{split=TRUE}.
#' 
#' @rdname markerPlots
#' @S3method profplot MarkerList
profplot.MarkerList <- function(x, y, groups=NULL, scale=FALSE, col=NULL, legend = FALSE, split = TRUE
								, ylab=NA, labels=NULL, ..., restore.gpar = TRUE){
	
	# extract expression values
	data <- y
	data <- exprs(data)
	# match markers
	i <- matchIndex(x, data)
	if( !nmark(i) ){
		warning("Nothing to plot: data does not contain any of the markers from the list.")
		return(invisible())
	}
	
	# drop empty indexes
	i <- drop(i)
	
	# change lengend specification into position if necessary
	if( !isFALSE(legend) ){
		if( isTRUE(legend) )
			legend <- 'topleft'
	}
	
	# define default ylab
	ylab_missing <- missing(ylab)
	
	if( isFALSE(split) ){ # all in one
		
		# extract types
		types <- as(dropvalues(i, int.ok=FALSE), 'factor')
#		str(types)
		
		# define colour if necessary
		if( is.null(col) )
			col <- rainbow(nlevels(types))
		col <- rep(col, length.out=nlevels(types))
		#
		# define default ylab
		if( ylab_missing ){
			ylab <- 'Marker expression values' 
			if( !isFALSE(scale) ) ylab <- paste(ylab, '(scaled)')
		}
		
		# do plot
		profplot(data[marknames(i), ], scale=scale, col=col[as.numeric(types)], ..., labels = labels, ylab=ylab, legend=FALSE)
		if( !isFALSE(legend) )
			legend(legend, legend=levels(types), col=col, lwd=1)
		
	}else{ # each type separately
		
		if( isTRUE(split) ) split <- mfrow(length(i))
		else if( !is.numeric(split) || length(split) != 2L )
			stop("Invalid argument `split`: should be a 2-long vector that indicates how to split plots (with mfrow).")
		
		if( length(split) > 1 ){
			op <- par(mfrow=split)
			if( restore.gpar ) on.exit( par(op) )
		}
		
		# setup scale
		if( isTRUE(scale) ) rep(TRUE, length(i))
		else if( isFALSE(scale) ) rep(FALSE, length(i))
		else if( is.numeric(scale) ){
			scale <- rep(scale, length.out=length(i))
			scale <- as.integer(pmin(scale, nmark(i, each=TRUE)))
		}else if( !is.list(scale) ) {
			stop('Invalid argument `scale`: must be either TRUE/FALSE, a numeric vector (integer) or a list.')
		}
		scale <- rep(scale, length.out=length(i))
		
		# extend legend
		legend <- rep(legend, length.out=length(i))
		
		res <- mapply(function(name, x, scale, legend, ...){
			# define colour if necessary
			if( is.null(col) ) col <- rainbow(length(x))
			col <- rep(col, length.out=length(x))
			#
			
			# subset the data
			data <- data[marknames(x),, drop=FALSE]

			# scale if necessary
			if( isTRUE(scale) ){
				data <- sweep(data, 2L, colMeans(data), '/')
				
			} else if( isInteger(scale) ){
				data <- sweep(data, 2L, data[scale,], '/')
				
			} else if( is.numeric(scale) ){
				scale <- rep(scale, length.out=ncol(data))
				data <- sweep(data, 2L, scale, '/')
			}
			# define default ylab
			if( ylab_missing ){
				ylab <- 'Expression values' 
				if( !isFALSE(scale) ) ylab <- paste(ylab, '(scaled)')
			}
			
			# do plot
			profplot(data, ..., col=col, legend=legend, main=name, ylab=ylab, labels=labels)
			
		}, names(i), i, scale, legend, MoreArgs=list(...))
		invisible(res)
	}
	
}

#' Cell Type Proportion Plot
#' 
#' Draws a vertical barplot of cell-types proportions.
#' 
#' @param x data.
#' Can be a matrix with cell types in rows and samples in columns, 
#' or an \code{\linkS4class{NMF}} object that contains proportions
#' as its mixture coefficient matrix, e.g. as computed by \code{\link{ged}}.
#' @param col color specification
#' @param legend.text logical that indicates if a legend should be added to 
#' the plot, or the text of the legend as in \code{\link{barplot}}.
#' @param scale logical that indicates if the data should be scaled so
#' that columns sum up to one, like proportions.
#' @param legwidth number of margin lines left to display the legend in the plot margin. 
#' @param ... extra arguments passed to \code{\link{barplot}}.
#' 
#' @export
propplot <- function(x, col=NULL, legend.text = TRUE, scale=FALSE, legwidth=10, ...){
	
	if( is.nmf(x) ) x <- coef(x)
	if( is.null(col) ) col <- rainbow(nrow(x))
	if( isTRUE(legend.text) ){
		legend.text <- paste(rownames(x)
				,'M', round(rowMeans(x, na.rm=TRUE), 2)
				,'S', round(apply(x, 1, sd), 2))
	}
	# plot
	if( scale ) x <- scoef(x)
	# setup layout for legend
	if( !is.null(legend.text) ){
		def.par <- par(no.readonly = TRUE)
		on.exit( par(def.par) )
		par(xpd=TRUE, mar=par()$mar+c(0,0,0,legwidth))
#		on.exit( close.screen(all.screens=TRUE) )
##		layout(matrix(c(1,2), nrow=1), width=c(0.9, .1))
#		split.screen(matrix(c(0,1,0,0.2,0,1,.2,1), nrow=2, byrow=TRUE))
#		screen(2)
	}
	res <- barplot(x, col = col, ...)
	if( !is.null(legend.text) ){
#		screen(1)
		legend('topleft', inset=c(1,0), legend=legend.text, fill=col)
	}
	# return res
	invisible(res)
}

emptyplot <- function(){
	plot(1, type = "n", axes=FALSE, xlab="", ylab="")
}


#' \code{screeplot} plots the following diagnostic plots:
#' \itemize{ 
#' \item for marker lists with associated scores, it plots the total number of 
#' markers whose values are over/under a range of threshold, as well as, 
#' optionally their distribtuion amongst the different cell types;
#' \item a screeplot of the condition number of a deconvolution matrix for a 
#' varying (increasing) number of markers.
#' }
#' 
#' @inheritParams stats::screeplot.default
#' @param breakdown logical that indicates if breakdowns of the number of markers
#' over the given range of values should be shown.
#' It may also be a single numberic that indicates the number of breakdowns to show.
#' @param range specifies the range of values to plot.
#' In \code{hist}, or in \code{screeplot} when \code{data=NULL}, it indicates a 
#' range over which the scores are used to computing the histogram. 
#' If \code{NULL}, then all values are used (Note: this might take a bit of time for long 
#' marker lists if \code{breakdown=TRUE}).
#' Otherwise it must be a 2-long numeric vector giving the interval of values to use.
#' 
#' If \code{data} is not \code{NULL}, it indicates which 
#' @param plot logical that indicates if the function should effectively produce a plot or 
#' only return the result.
#' 
#' Otherwise, if \var{data} is not \code{NULL}, then it must be either a range of values 
#' 
#' @return When \code{data=NULL} and \code{breakdown} is not \code{FALSE} this function returns a 
#' list with the result from the breakdown barplot and the breakdown (in element \sQuote{breakdown}).
#' 
#' When \code{data} is not \code{NULL}, a \code{\linkS4class{MarkerList}} object made of the set of markers
#' that achieve the least condition number is returned.
#' 
#' @rdname markerPlots
#' @S3method screeplot MarkerList 
screeplot.MarkerList <- function(x, data=NULL, breakdown=FALSE, range=NULL
								, xlab=NULL, ylab=NULL, main=NULL, ..., plot=TRUE){
	
	if( !is.null(data) ){
		
		nm <- nmark(x)
		if( !nm ){
			warning("Nothing to plot: the marker list object is empty.")
			return(invisible())
		}
		
		r <- kappa(x, data, n=range)
		# extract minimum
		mins <- attr(r, 'mins')
		nopt <- attr(r, 'n')[1L] 
		s <- flatten_seq(x, n=nopt)
		
		if( missing(xlab) ) xlab <- 'Number of markers'
		if( missing(ylab) ) ylab <- 'Condition number'
		if( missing(main) ) main <- 'Signature matrix condition number'

		# plot if requested
		if( plot ){
			plot(names(r), r, ..., xlab=xlab, ylab=ylab, main=main)
			abline(v=nopt, col='red')
			axis(1, at=nopt, labels=nopt)
		}
		res <- attr(s, 'markers')
		attr(res, 'kappa') <- r
		# return result
		if( plot ) return(invisible(res))
		else return(res)
		
	}
	
	# check for values
	if( !hasValues(x) ){
		warning("Nothing to plot: the marker list object does not contain any numeric values/scores.")
		return(invisible())
	}
	
	# compute range 
	v <- unlist2(x)
	
	if( !is.null(range) && (!is.numeric(range) || length(range) != 2L) ){
		stop("Invalid range: if not NULL, should be of length 2 = c(min, max) [", class(range), "]")
	}
	if( !isFALSE(breakdown) ){
		
		if( isTRUE(breakdown) ) breakdown <- min(nmark(x)/10, 20)
		if( is.null(range) ) range <- c(min(v), max(v)) 
		range <- seq(range[1], range[2], length.out = breakdown)
		
		# compute breakdown of markers for each value in range
		sdata <- sapply(range, function(lambda){
			v <- v[v<=lambda]
			x <- x[,names(v)]
			sapply(x, length)
		})
	
		# plot
		tot <- colSums(sdata)
		r <- round.pretty(range, 2)
		b <- barplot(sdata, names.arg=paste(r, '\n(', tot, ')', sep=''), beside=TRUE, ..., plot=plot)
		
		# combine result from barplot and breakdown
		res <- c(b, list(breakdown=sdata))
		# return result
		if( plot ) return(invisible(res))
		else return(res)
		
	}else{
		# plot histogram
		hist(x, range=range, ...)
	}
}

#' \code{hist} plots the histogram of numeric scores associated with each marker in a marker list, 
#' if any.
#' 
#' @inheritParams graphics::hist
#' @S3method hist MarkerList
#' @rdname markerPlots
hist.MarkerList <- function(x, range=NULL, split=FALSE, xlab='Marker scores', main='Histogram of marker scores', ..., restore.gpar = TRUE){
	
	# check for values
	if( !hasValues(x) ){
		warning("Nothing to plot: the marker list object does not contain any numeric values/scores.")
		return(invisible())
	}
		
	if( !is.null(range) && (!is.numeric(range) || length(range) != 2L) ){
		stop("Invalid range: if not NULL, it should be of length 2 = c(min, max) [", class(range), "]")
	}
	
	if( isFALSE(split) ){ # all in one
		# unlist
		v <- unlist2(x)
		
		# limit to range if necessary
		if( !is.null(range) ){
			v <- v[ v >= range[1] & v <= range[2]]
		}
		
		# plot
		hist(v, ..., xlab=xlab, main=main)
		
	}else{ # each type separately
		
		# remove empty sets
		x <- drop(x)
		
		# determine how to split 
		if( isTRUE(split) ) split <- mfrow(length(x))
		else if( !is.numeric(split) || length(split) != 2L )
			stop("Invalid argument `split`: should be a 2-long vector that indicates how to split plots (with mfrow).")
		
		if( length(split) > 1 ){
			op <- par(mfrow=split)
			if( restore.gpar ) on.exit( par(op) )
		}
		
		res <- sapply(names(x), function(i, ...){ hist(x[i], range=range, main=i, ..., xlab=xlab) }, ...)
		invisible(res)
	}
}

round.pretty <- function(x, min=2){
	
	if( is.null(x) ) return(NULL)		
	n <- 0
	y <- round(sort(x), n)
	while( any(diff(y)==0) ){
		n <- n+1
		y <- round(sort(x), n)
	}	
	round(x, max(min,n))
}


#' Condition Number of a Marker List
#' 
#' The S3 method \code{kappa} for \code{MarkerList} objects returns the 
#' condition number of a data (sub-)matrix limited to markers in the list.
#' 
#' This method typically is used to predict/optimise the deconvolution 
#' power of a set of cell-specific signatures, as proposed by \cite{Abbas2009}.
#' 
#' @inheritParams base::kappa
#' @param z a \code{\linkS4class{MarkerList}} object
#' @param data reference data matrix or \code{ExpressionSet} object, 
#' from which is extracted the sub-matrix of markers.
#' @param n specifies the range of total number of markers that compose  
#' the matrix whose condition number is computed:
#' \itemize{
#' \item The default (\code{0}) is to compute using all makers in \var{z} that 
#' can be found in the data.
#' \item A single positive value indicates the exact number of such markers to use.
#' A warning is thrown if this number is greater than the number of markers
#' present in the data. 
#' \item A single negative value does the same as a positive value, but no warning is 
#' thrown if it exceeds the actual number of markers in the data.
#' \item \code{NULL} indicates a full range of number of markers, i.e. the computation is 
#' performed including for matrices composed of \code{length(z)} to \var{n} markers.
#' \item A numeric vector completely specifies the range of number of markers over which to 
#' perform the computation.
#' }
#' 
#' @param ... extra arguments passed to \code{\link{kappa}}.  
#' 
#' @seealso \code{\link{kappa}}
#' @S3method kappa MarkerList
#' @examples 
#' 
#' # random data and markers
#' d <- rmatrix(20, 10)
#' x <- rMarkerList(4, 5)
#' summary(x)
#' kappa(x, d)
#' 
#' # over a specifc or range of total number of markers
#' kappa(x, d, 6)
#' kappa(x, d, 1:20)
#' 
#' # condition number is Inf if no markers are found
#' # NB: throws a warning
#' rownames(d) <- as.character(1:nrow(d))
#' x <- rMarkerList(4, 5, names=TRUE)
#' kappa(x, d)
#' 
#' 
kappa.MarkerList <- function(z, data, n=0, ...){
	
	# drop emtpy types
	x <- drop(z)
	# match index
	i <- matchIndex(x, data)
	# early exit if no markers matched
	if( !nmark(i) ){
		warning("Could not find any markers in data")
		return( Inf )
	}
	nm <- nmark(i)
	if( isNumber(n) ){
		if( n < abs(n) ){
			if( n > 0 )
				warning("Using all ", nm," markers to compute the condition number (asked for ", abs(n), ").")
			n <- 0
		}
		n <- abs(n)
	}
	
	if( is.null(n) ) n <- 1:nm
			
	# extract expression data
	data <- exprs(data)
	if( identical(n, 0) || identical(n, 0L) ){
		kappa(data[marknames(i), ], ...)
	}else{
	
		# re-shape the list to simplify its parsing
		idx <- flatten_seq(i, decreasing=FALSE)
		if( isNumber(n) ){
			return( setNames(kappa(data[head(idx, n), ], ...), n) )
		}
		
		# limit to relevant range value
		n <- n[n>=length(x)]
		
		# initialise sub-matrix
		r <- sapply(n, function(r, ...) { kappa(data[head(idx, r), ], ...) }, ...)
		res <- setNames(r, n)
		
		# add attributes for minimum
		imin <- which(res == min(res))
		mins <- res[imin]
		attr(res, 'mins') <- mins
		attr(res, 'n') <- as.numeric(names(mins))
	
		res
	}
	
	
}

#' Heatmaps Highlighting Markers
#' 
#' The function \code{markermap} draws a heatmap of a reference expression data 
#' (e.g., the expression matrix from pure samples or estimated cell-specific signatures), 
#' where marker are annotated by colored bands on the left-hand side of the heatmap.
#' 
#' Argument \code{view} controls the way markers are annotated.
#' On all views, markers are coloured according to their type, defined by the 
#' element of \code{object} in which they appear.
#' Each couloured tick/cell corresponds to a different marker position.
#' 
#' The following views are available:
#' \describe{
#' \item{single}{ a single row annotation is added, showing the position of
#' each marker.}
#' \item{split}{ one track per marker type in added, showing the position of 
#' each marker in its respective cell type.}
#' \item{predict}{ one track per column in \code{data} is added, showing the 
#' position of each marker in the most expressing column.
#' When \code{data} is a basis matrix obtained from deconvolution, this view 
#' is useful to check how known cell type markers (the coulours) map on 
#' estimated signatures (the annotation columns).}
#' } 
#' 
#' @param object a \code{MarkerList} object
#' @param data reference data object, whose values are used in the heatmap.
#' @inheritParams NMF::aheatmap
#' @inheritParams .atrack,MarkerList-method
#' @param subsetRow this argument acts as in \code{\link{aheatmap}}, but 
#' if \code{subsetRow=TRUE}, then the heatmap is limited to the markers
#' only.
#' 
#' @inline
#' @export
#' 
#' @examples 
#' 
#' x <- rmix(3, 100, 20)
#' m <- getMarkers(x)
#' markermap(m, basis(x))
#' markermap(m, x, view='single')
#' basismarkermap(m, rnmf(3, x))
#' 
#' # after real deconvolution 
#' res <- ged(x, coef(x), 'csSAM')
#' basismarkermap(m, res)
#' markermap(m, res, view='split')
#' 
setGeneric('markermap', function(object, data, ...) standardGeneric('markermap') )
#' Workhorse method is for \code{markermap}.
setMethod('markermap', signature(object='MarkerList', data='matrix'), 
	function(object, data, annCol=NA, annColors=NA, annRow=NA, view='split', subsetRow=NULL, scale='row' 
				, color='YlOrRd:100', Rowv=TRUE, Colv=NA, distfun="correlation", hclustfun="average"
				, ...){
		
		# limit to markers
		if( isTRUE(subsetRow) || identical(subsetRow, 'markers') ){
			i <- matchIndex(object, data)
			data <- data[marknames(i), ]
			subsetRow <- NULL
			# remap if integer marker list
			if( mltype(object, 'integer') ){
				j <- 0L
				object <- sapply(object, function(x){
							i <- seq_along(x) + j
							j <<- j + length(x)
							setNames(i, names(x))
						}
					, simplify=FALSE)
			}
		}
		#
		
		# add matching column annotation  
		if( isTRUE(annCol) ){
			if( setequal(colnames(data), names(object)) ) annCol <- colnames(data)
			else if( is.null(colnames(data)) && length(object) == ncol(data) ) annCol <- names(object)
			else annCol <- NA
		}
		
		# construct feature annontation data.frame
		rann <- .atrack(object, data=data, view=view)
		lev <- unique(unlist(lapply(rann, levels)))
		lev <- lev[match(names(object), lev)]
		annRow <- atrack(annRow, rann, .DATA=NMF:::amargin(data, 1L))
		#str(lev)
		annColors <- c(annColors, list(Markers=setNames(rainbow(length(lev)), lev)))
		
		# generate the heatmap
		aheatmap(data, scale=scale, Rowv=Rowv, Colv=Colv, distfun=distfun, hclustfun=hclustfun 
				, subsetRow = subsetRow, annCol=annCol, annRow=annRow
				, annColors=annColors
				, color=color, ...)
	}
)

#' The method \code{markermap} for \code{ExpressionSet} objects calls the main 
#' \code{markermap} method on the expression matrix \code{exprs(object)}.
#' 
setMethod('markermap', signature(object='MarkerList', data='ExpressionSet')
		, function(object, data, ...){
			markermap(object, exprs(data), ...)
		}
)
#' The method \code{markermap} for \code{NMF} objects calls the main 
#' \code{markermap} method on the basis matrix \code{basis(object)}. 
#' 
setMethod('markermap', signature(object='MarkerList', data='NMF')
	, function(object, data, ...){
		markermap(object, basis(data), ...)
	}
)
#' The method \code{markermap} for \code{NMFfitX} objects calls the main 
#' \code{markermap} method on the best fit \code{fit(object)}. 
#' 
setMethod('markermap', signature(object='MarkerList', data='NMFfitX')
	, function(object, data, ...){
		
		markermap(object, fit(data), ...)
	}
)
#' This method extracts and plots a list of markers from 
#' a set of basis signature matrix (in \code{object}), 
#' where each feature is associated with the most-expressing signature.
#' 
setMethod('markermap', signature(object='MatrixData')
	, function(object, data, ...){
		ml <- MarkerList(.extractMarkers_signatures(object))
		markermap(ml, data, ...)
	}
)

#' \code{basismarkermap} calls \code{markermap} with arguments tuned so that: 
#' no column reordering is performed, the rows are scaled to sum up to one, the heatmap 
#' only shows marker expression values, and the markers are placed according to their 
#' most expressing estimated signature.
#' This view is meant to help in assessing the validity of deconvolution results, when 
#' known markers are available: ideally marker row annotations should be composed of monochrome 
#' monoblocks.
#'  
#' @rdname markermap
#' @export
basismarkermap <- function(object, data, scale='r1', view='predict', subsetRow=TRUE, Rowv=NA, labRow=NA, ...){
	
	markermap(object, data, view=view, ..., subsetRow=subsetRow, Rowv=Rowv, labRow=labRow, scale=scale)
}


#' Plotting MarkerList Objects
#' 
#' The function \code{barplot} plots a barplot of the number of markers of each cell type
#' present in a \code{\linkS4class{MarkerList}} object, or their respective expression values 
#' in a given dataset.
#' 
#' @inheritParams graphics::barplot
#' @param height a \code{\linkS4class{MarkerList}} object.
#' @param data dataset (a matrix or \code{\link{ExpressionSet}} object) that contains expression 
#' values of the markers in some experiment.
#' If not missing, values corresponding to markers are plotted for each column in the data. 
#' Markers that do not have corresponding values in \code{data} are discarded.
#' Note that no identifier conversion is performed, meaning that the marker list and the dataset 
#' must use the type of gene identifiers.
#' Use the function \code{\link{convertIDs}}, if a conversion is needed.
#' @param col colour specification for the bars.
#' @param byType logical that indicates if the colour specification in \code{col} should be 
#' interpreted as colours for each cell type.
#' Is is used only when \code{data} is not missing.
#' If \code{TRUE} then \code{col} specifies the colours to be used for all bars in each cell type.
#' If \code{FALSE} 
#' @param ylab label for the y-axis. Default labels are provided if missing.
#' @param scale logical that indicates if the data should be scaled. 
#' This argument is used only when \code{beside=FALSE}.
#' If \code{TRUE} the gene expression profiles (rows of \code{data}) are scaled to sum up to 100,
#' then the sample expression profiles (the columns of \code{data}) are scaled to sum up to 100.
#' The scaled values can be interpreted as how much a marker gene is expressed in a sample 
#' relatively to other samples and genes.
#' @param legend logical that indicates if a legend box should be added.
#' @param reorder when \code{data} is missing, this argument indicates if/how the cell types 
#' should be reordered: \code{FALSE} do not reorder the cell types, a negative number orders 
#' by decreasing number of markers, and anything other than \code{FALSE} orders by increasing order 
#' of number of markers.
#' 
#' When \code{data} is not missing, then it indicates if the markers should be ordered in the same order as 
#' they appear in the data. 
#' 
#' is used when \code{data} is not missing and indicates if 
#' the markers should be ordered as they appear in the dataset, as opposed to keep their 
#' original order, i.e. grouped by cell type.
#' @param las specifies the orientation of the x-labels: 1 for horizontal, 2 for vertical. 
#' @param ... extra arguments passed to \code{\link[graphics]{barplot}}.
#' 
#' @importFrom graphics barplot
#' @S3method barplot MarkerList
#' 
#' @export
barplot.MarkerList <- function(height, data=NULL
								, col=NULL, byType=TRUE
								, ylab=NULL, border=missing(data)
								, beside=TRUE, scale=!beside, legend=TRUE
								, reorder=if( is.null(data) ) -1L else FALSE, las=2, ...){
		
		x <- height
		y <- data
		n <- nmark(x, each=TRUE)
		n.min <- min(n)
		n.max <- max(n)
				
		# draw barplot of number of markers if no data is provided
		if( is.null(y) ){
			nm <- nmark(x, each=TRUE)
			res <- list()
			idx <- NULL
			
			if( !isFALSE(reorder) ){
				dec <- if( isNumber(reorder) ) reorder<0 else TRUE
				if( !is.null(dec) ){
					idx <- order(nm, decreasing=dec)
					nm  <- nm[idx]
					# reorder the colors
					if( length(col) == length(nm) ) col <- col[idx]
				}
			}
			# add default label to y-axis
			if( missing(ylab) )
				ylab <- "Number of markers"
			# draw the plot
			res$midpoint <- barplot(nm, beside=beside, ylab=ylab, col=col, las=las, ...)
			res$ylim <- c(0, max(nm))
			res$idx <- idx
			res$col <- col
			return(invisible(res))
			
		}else{
			
			# extract indices
			lidx <- matchIndex(x, y)
			nm <- nmark(lidx, each=TRUE)
			lidx.ok <- nm>0L
			# drop empty types
			lidx <- drop(lidx)
			if( !length(lidx) ){
				warning("barplot - None of the markers were found in the data: nothing to plot.", call.=FALSE)
				return(invisible())
			}
			if( isTRUE(byType) ){ # interpret colours as one for each cell-type
				
				if( is.null(col) ){
					col <- rainbow(length(lidx))
				}else{
					# convert col if nessecary 
					if( is.logical(col) ) col <- as.factor(col)
					if( is.factor(col) ) col <- grey.colors(nlevels(col))[as.integer(col)]
					
					if( length(col) < length(x) ){
						# expand `col` if necessary 
						col <- rep(col, length.out=length(x))
					}else if( length(col) > length(x) ){
						warning("Argument `col` has more elements [", length(col) ,"] than `x` [", length(x),"]:"
									, " discarding ", length(col) > length(x), " colour(s).")
						col <- head(col, length(x))
					}
				}
				stopifnot( length(col) == length(x) )
				# subset the colours
				if( length(col) == length(x) ) col <- col[lidx.ok]
				# store type colors 
				tcol <- col
				
				col <- unlist(mapply(rep, col, nmark(lidx, each=TRUE), SIMPLIFY=FALSE))
				
			}else if( !is.null(col) ){
				col <- rep(col, length.out=nmark(lidx))
			}
			
			# extract data
			edata <- if( is(y, 'ExpressionSet') ) exprs(y) else as.matrix(y)
			
			idx <- unlist2(lidx)			
			if( !is.null(col) && length(col) != length(idx) ){
				stop("Number of colours [", length(col) ,"] is not equal to the number of markers [", length(idx),"]")
			}
			
			# if requested: reorder in the same order as in the expression data
			# as opposed to the order within each cell-type
			if( reorder ){
				o_idx <- order(idx)
				if( length(col) == length(idx) ) col <- col[o_idx]
				idx <- idx[o_idx]
			}
			
			# create subset data
			.data <- edata[idx,,drop=FALSE]
			
			if( !beside && scale ){
				#print(length(unlist(mapply(rep, nm, nm, SIMPLIFY=FALSE))))
				#.data <- sweep(.data, 1L, unlist(mapply(rep, nm, nm, SIMPLIFY=FALSE)), '/')
				.data <- sweep(.data, 1L, rowSums(.data), '/') * 100
				.data <- sweep(.data, 2L, colSums(.data), '/') * 100
			}
			
			# add default label to y-axis
			if( missing(ylab) )
				ylab <- "Marker expression value"
			
			
			# plot barplot
			mp <- barplot(.data, col=col, beside=beside, ylab=ylab, border=border, las=las, ...)
			if( !isTRUE(byType) )
				return(invisible(list(midpoint=mp, ylim=c(0, max(.data)), idx=lidx)))
			
			# draw legend if requested
			legpos <- 'topleft'
			if( is.character(legend) ){
				legpos <- legend
				legend <- TRUE
			}
			
			if( !isFALSE(legend) ){
				nm.ok <- nmark(lidx, TRUE)
				l <- length(lidx)
				legend(legpos, legend=paste(names(lidx), ' (', nm.ok, ')', sep='')
					, fill=tcol
					, ncol= if(l > 4) ceiling(l/4) else 1
					, inset=0.05)
			}
			# return the bars' mid-positions from barplot 
			invisible(list(midpoint=mp, col=tcol, idx=which(lidx.ok), ylim=c(0, max(.data))))
			
		}

		
	}
#)


#' Split Boxplot by Group
#' 
#' Plots a boxplot with each data split into sub-boxplots, defined 
#' according a factor variable.
#' 
#' @param x data matrix, with samples in column.
#' 
#' NB: in \code{\link{boxplot}} \var{x} has samples in rows.
#' 
#' @rdname boxplotBy
#' @export
boxplotBy <- function(x, ...){
	UseMethod('boxplotBy')
}
#' The default method plots each \strong{column} using one boxplot for each 
#' level in a factor.
#' 
#' @param by grouping variable, coerced to a factor if necessary.
#' It is recycled if not of length the number of columns in \var{x} 
#' @param col spcecifies the colour to use for each group by argument \var{by}.
#' @param legend legend specification.
#' If a logical, it turns on/off the legend, if a character string, it indicates
#' the position of the legend, e.g. \code{'topright'} (see \code{\link{legend}} 
#' for possible values).
#' @inheritParams graphics::boxplot
#' @param ... extra arguments passed to the base \code{\link{boxplot}} function.
#'   
#' @author This function was adapted from the following post on R-help:
#' 
#' \url{https://stat.ethz.ch/pipermail/r-help/2003-March/031046.html}
#' 
#' @S3method boxplotBy default
#' @rdname boxplotBy
#' @examples 
#' 
#' x <- rmatrix(3, 20)
#' boxplotBy(x, gl(2, 10))
#' 
boxplotBy.default <- function(x, by, scale=FALSE, col=rainbow(nlevels(by)), legend=TRUE, names=colnames(x), ...){
		
	# convert groups into a factor
	if( !is.factor(by) ) by <- factor(by, levels=unique(by))
	groups <- by
	
	df <- as.data.frame(x)
	# extend groups
	nvar <- ncol(df)
	groups <- rep(groups, length.out=nrow(df))
	
#	df <- data.frame(v1 = rnorm(100), v2 = rnorm(100), 
#			v3 = rnorm(100), grp = rep(letters[1:3],100))
#	
	
# Define 'at' argument positions. There will be 9 boxes generated
# by the boxplot() call normally centered on x axis1:9. 
# Using 'at', you can reconfigure these.
# Note that the center boxes for each variable
# are at x axis 2, 5 and 8. The other two bars are offset by 0.5
# before and after the center box.
	
	n <- nlevels(groups)
#	at = c(1.5, 2.0, 2.5, 4.5, 5.0, 5.5, 7.5, 8.0, 8.5)
	at <- seq(1, nvar*n)/n
	shift <- unlist(mapply(rep, seq(1,nvar), rep(n, nvar), SIMPLIFY=FALSE))
	at <- at + shift
# Generate plots as per Thomas' example, cycling colors and box
# positions. Set boxwex for thin boxes. Do not draw axes.
	bp <- boxplot(do.call("c",lapply(df,function(v) split(v,groups))), 
			col = rep(col, length.out=length(at)), boxwex = 0.2, at = at, 
			axes = FALSE, ...)
	
# Draw y axis
	axis(2)
	
# Draw x axis, creating letter labels below each box
	#axis(1, labels = rep(levels(groups), ncol(df)), at = at, las=2)
	
	# legend if requested
	legpos <- 'topleft'
	if( isString(legend) ){
		legpos <- legend
		legend <- TRUE
	}
	if( isTRUE(legend) )
		legend(legpos, fill=col, legend=levels(groups))
	
# Now label each group of 3 boxes with the varname in the middle
# of each grouping.
	
	mid <- as.numeric(by(at, shift, mean))	
	axis(side = 1, names, at = mid, las=2)
	
# Draw box around the whole plot
	
	box()
	
	# return plot data
	invisible(bp)
		
}

#' The method for \code{\linkS4class{NMF}} objects plots use the data from the 
#' coefficient matrix. 
#' 
#' @param scale logical that indicate if each data column should be scaled so 
#' that it sums up to one.
#' @rdname boxplotBy
#' @S3method boxplotBy NMF
boxplotBy.NMF <- function(x, by, scale=FALSE, ...){
	
	# extract coef
	x <- coef(x)
	# scale coef 
	if( scale ){
		x <- pmax(x, 0)
		x <- sweep(x, 2L, colSums(x), '/')
	}
	boxplotBy.default(x=t(x), by=by, ...)
}

#' Plots Cell-Specific FDR Estimates
#' 
#' @param x data object, typically returned by \code{\link{ged}}.
#' @param ... extra parameters passed to subsequent calls.
#' @param filename file where to save the plot.
#' Can be of any type amongst pdf, svg, png, bmp, jpg, tiff.
#' 
#' @export
csplot <- function(x, ..., filename=NULL){
	if( !is.null(filename) ){
		NMF:::gfile(filename)
		on.exit( dev.off() )
	}
	UseMethod('csplot')
}

#' @S3method csplot NMFfit
#' @rdname csplot
csplot.NMFfit <- function(x, ...){
	x <- structure(list(x), class=str_c('ged_', algorithm(x)))
	csplot(x, ...)
}

#' @S3method csplot NMFfitX
#' @rdname csplot
csplot.NMFfitX <- function(x, ...){
	csplot(minfit(x), ...)
}
