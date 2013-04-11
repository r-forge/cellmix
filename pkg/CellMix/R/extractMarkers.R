# TODO: Add comment
# 
# Author: renaud
# Created: 28 Nov 2012
###############################################################################

#' @include markers.R
NULL

#' Extract Markers from Pure Samples.
#' 
#' When expression data from pure samples are available, it is possible to 
#' select genes for their ability to discriminate between the different cell-types.
#' The function \code{extractMarkers} provides an interface to compute scores or statistics 
#' using different methods, and select the ones that pass a given threshold.
#' 
#' To see all available methods use \code{markerScoreMethod()}.
#' 
#' @inline
#' @export
#' @examples
#' # see available scoring methods
#' markerScoreMethod()
#' 
#' 
#' 
setGeneric('extractMarkers', function(object, ...) standardGeneric('extractMarkers') )

#' @param object A numeric \code{matrix}, an object of class 
#' \code{\link[Biobase]{ExpressionSet}}, or a \code{\linkS4class{MarkerList}} 
#' object. 
#' 
#' @param data main extra data used by the scoring and selection methods.
#' 
#' If \var{object} is a matrix-like object, \var{data} is generally a factor or 
#' grouping variable, as a vector, that defines the 
#' cell-type for each -- pure -- sample. If a vector, this argument is 
#' converted into a factor with the levels in their order of appearance in 
#' \code{data}, by \code{factor(data, levels=unique(data))}. 
#' This is to obtain levels in an order that is consistent with the samples' order.
#' 
#' If \var{object} is a \code{MarkerList} object, then \var{data} is generally a 
#' matrix-like object that contains expression data.
#' 
#' @param method accession key of the method either used to compute marker scores 
#' in \code{extractMarkers}, or to define/retrieve scoring methods in \code{markerScoreMethod}.
#' 
#' All available methods can be retrieved via \code{markerScoreMethod()}.
#' 
#' A custom scoring method can also be passed as a function.
#' @inheritParams selectMarkers
#' @param ... other parameters passed to the scoring method in \code{extractMarkers} \strong{and} \code{\link{selectMarkers}}, 
#' or that define default arguments when defining a scoring method with \code{markerScoreMethod}.
#' @param format specifies the output format: \sQuote{list} returns a 
#' \code{\linkS4class{MarkerList}} object and  \sQuote{raw} directly returns the 
#' result of the scoring method -- which can be useful for lengthy computations.
#'  
#' 
setMethod('extractMarkers', 'ANY', 
		function(object, data=NULL, method = 'Abbas', threshold=NA, decreasing=FALSE
				, ..., format=c('list', 'raw')){
			
			# check validity of argument 'data'
			if( is.matrix(object) && !is.null(data) && (is.vector(data) || is.factor(data)) ){
				if( length(data) != ncol(object) )
					stop("extractMarkers - Incompatible arguments: the length of `data` must be equal to the number of columns of `object`.")
				
				# convert 'data' into a factor if necessary
				if( !is.factor(data) ) data <- factor(data, levels=unique(data))
				data <- droplevels(data)
			}
			
			# match format
			format <- match.arg(format)
			
			# get method from registry
			f <- markerScoreMethod(method, class(object))
			if( is.list(f) ){
				# load defaults from method
				method <- f$key
				# load default argument values
				if( missing(decreasing) && !is.null(f$decreasing) ) decreasing <- f$decreasing
#				if( missing(statistic) && !is.null(f$statistic) ) statistic <- f$statistic
				if( missing(threshold) && !is.null(f$threshold) ) threshold <- f$threshold
				# extract function (first element)
				f <- f[[1L]]
				if( !is.function(f) )
					stop("Invalid registered marker scoring method: first element should be a function [", class(f), ']')
			}
			if( !is.function(f) )
				stop("Invalid marker scoring method: must be a function [", class(f), ']')
			
			useData <- hasFormals(f, 'data') && (!is.null(data) || hasFormals(f, 'data', withdefault=FALSE) ) 
			# compute scores
			if( useData ){
				stats <- f(object, data=data, ...)
				# add the grouping variable as an attribute
				attr(stats, 'types') <- levels(data)
			}else{
				if( !is.null(data) ){
					warning("Discarding argument `data`: scoring method "
							, if( isString(method) ) str_c(sQuote(method), ' ')
							, "has no `data` argument.\n"
							, "  Method signature: ", str_fun(f))
				}
				stats <- f(object, ...)
			}
			
			if( !isS4(stats) ){
				cl <- 'markerScore'
				if( isString(method) ) cl <- c(str_c('markerScore_', method), cl)
				class(stats) <- c(cl, class(stats))
			}
			
			# return plain result if requested
			if( format == 'raw' ) return(stats)
			
			# select markers based on scores
			selfun <- selectS3method('selectMarkers', class(stats), optional=FALSE)
			if( is.null(selfun) ){
				stop("Could not find a suitable marker selection method for result object of class(es) ", str_out(class(stats), Inf))
			}
			ca <- call('selfun', stats, threshold = threshold, decreasing = decreasing, quote(...))
#			if( hasFormals(selfun, 'statistic') ) ca$statistic <- quote(statistic)
			if( useData && hasFormals(selfun, 'data') ) ca$data <- quote(data)
			if( hasFormals(selfun, '.object') ) ca$.object <- quote(object)
			res <- eval(ca)
			
			return( res )
		} 
)

#' The method \code{extractMarkers} for \code{ExpressionSet} objects calls 
#' \code{extractMarkers} on the expression matrix \code{exprs(object)}.
#' 
setMethod('extractMarkers', 'ExpressionSet'
	, function(object, ...){
		res <- extractMarkers(exprs(object), ...)
		# pass on annotations
		if( isMarkerList(res) && hasAnnotation(object) ){
		  annotation(res) <- annotation(object)
		}
		res
	}
)

#' Select Markers Based on Scores 
#' 
#' The function \code{selectMarkers} filters raw results returned by scoring methods  
#' and wrap them into a \code{MarkerList} object.
#' It is called by \code{\link{extractMarkers}} when \var{format} is not \code{'raw'}.
#' 
#' @inheritParams extractMarkers
#' @param x data object based on which the markers are selected, 
#' as computed returned by \code{extractMarkers(..., format='raw')}.
#' The type of \code{x} depends on the scoring method used to compute it.
#' @param .object argument used internally by \code{\link{extractMarkers}} 
#' to pass the original object on which marker scores where computed.
#' 
#' @return \code{selectMarkers} returns an object of class \code{MarkerList}.
#' @export
selectMarkers <- function(x, ..., .object=NULL){
	UseMethod('selectMarkers')
}

#' @param threshold threshold that is applied to filter markers based on their respective 
#' statistic/score: 
#' Genes with \code{statistic >= threshold} are selected if \code{decreasing=TRUE}.
#' Otherwise, if \code{decreasing=FALSE} or \code{NA}, the selected genes are those with 
#' \code{statistic <= threshold}.
#' 
#' Filtering is disabled if \code{threshold=NA}.
#' @param decreasing logical that indicates how the statistic/score should be 
#' ordered: \code{TRUE} orders by decreasing value (i.e. the greater the score the better),
#' while \code{FALSE} orders by increasing value (i.e. the lower the score the better).
#' The value of this argument also affects the way the filtering is performed (see argument 
#' \code{threshold}).
#' 
#' Ordering is disabled if \code{decreasing=NA}.
#' 
#' @rdname selectMarkers
#' @S3method selectMarkers MarkerList
selectMarkers.MarkerList <- function(x, threshold=NA, decreasing=FALSE, ..., .object){
	
	# apply threshold if necessary
	if( !is_NA(threshold) ){
		sel <- if( is_NA(decreasing) || decreasing ){
			x >= threshold
		}else{
			x <= threshold
		}
		x <- x[ sel ]
	}
	# reorder if necessary
	if( !is_NA(decreasing) ){
		x <- reorder(x, decreasing = decreasing)
	}
	# return filtered/reordered object
	x
}
 
#' @rdname selectMarkers
#' 
#' @details
#' The default method \code{selectMarkers.markerScore} selects markers based on 
#' scores stored as a matrix-like object, whose first column contains a factor or 
#' a vector that assigns each feature to a group (cell type) and column with index 
#' \code{statistic} (by default the second column) contains the feature scores.
#' 
#' Scoring methods that do not define a dedicated \code{selectMarkers} must return a 
#' \code{MarkerList} object or a matrix that complies with the description above. 
#' It is recommended that, if such methods compute a range of statistics, they 
#' store the requested selection statistic in the second column.
#' 
#' @param statistic name or column index of the statistic/score to use.
#' The default is to use the second column of the matrix returned by the scoring
#' method, but each method can define its own default.
#' 
#' @S3method selectMarkers markerScore
selectMarkers.markerScore <- function(x, statistic=2L, data=attr(x, 'types'), ..., .object){
	
	stats <- x
	itype <- stats[, 1L]
	if( is.null(data) ){
		if( is.numeric(itype) )
			stop("missing argument `data`: argument required when the input data has no attribute `types`.")
		
		# first column contains the types
		group <- itype
		if( !is.factor(group) )
			group <- factor(itype, levels=unique(itype))
		# wrap into a MarkerList object
		m <- setNames(stats[, statistic], rownames(stats))
		m <- MarkerList(split(m, group))
	}else{
		
		# extract levels if types is a factor
		if( is.factor(data) ){
			data <- droplevels(data)
			data <- levels(data)
		}
		
# 		g <- unique(itype)
# 		if( length(data) != length(g) )
# 			stop("number of types in `data` [", length(data), "] must be equal to the number of groups of scores [", length(g), ']')
		
		# wrap into a MarkerList object
		m <- setNames(stats[, statistic], rownames(stats))
		m <- MarkerList(split(m, data[itype]))
	}
	
	selectMarkers(m, ...)
}

#############################################################
#
#############################################################
# register using pkgmaker registry features
GEDscore_registry <- registry::registry(registry_class = "GEDscore_registry"
		, entry_class = "GEDscore_entry")

# access key
GEDscore_registry$set_field("key", type="character"
		, is_key = TRUE
		, index_FUN = match_partial_ignorecase
		, validity_FUN = checkKey)
# type of data object
GEDscore_registry$set_field("object", type="character", mandatory = TRUE)
# method
GEDscore_registry$set_field("algorithm", type="list", mandatory = TRUE)

GEDscore_registry <- setPackageRegistry('score', GEDscore_registry
		, description = 'Methods to score and select marker genes'
		, entrydesc = 'marker scoring method')

#' \code{markerScore} is a shortcut for \code{extractMarkers(..., format='raw')}.
#' 
#' For long-runnning scoring methods, the marker selection process could be performed manually, 
#' by separating the scoring and selections steps, e.g.:
#' \samp{
#' sc <- scoreMarkers(eset, data=groups, method='HSD')
#' selectMarkers(sc, threshold=10^-3)
#' }
#'   
#' @family markerScore
#' @rdname extractMarkers 
#' @export
scoreMarkers <- function(...) extractMarkers(..., format='raw')

#' \code{markerScoreMethod} defines/retrieves marker scoring methods.
#' 
#' @family markerScore
#' @rdname extractMarkers 
#' @export
markerScoreMethod <- function(method, object, ...){
	if( missing(method) ) GEDscore_registry$get_entry_names()
	else if( nargs() <= 2L ){
		if( isString(method) ){
			obj <- pkgreg_fetch('score', key=method)
			if( !object %in% obj$object ){
				stop("Scoring method ", sQuote(method), " cannot handle objects of type '", object, "'"
					, " [types: ", str_out(obj$object, Inf),"]")
			}
			c(obj$algorithm, key=obj$key)
		}else if( is.function(method) ) method
		else if( is.list(method) && length(method) && is.function(method[[1L]]) ) method
		else stop("Scoring method must be specified as a string, a function or a list whose first element is a function [", class(method), "]")
	}else{
		args <- list(...)
		if( !identical(args, list(NULL)) ){
			
			ifun <- 1L
			if( !is.null(na <- names(args)) ) 
				ifun <- which(na=='')[1L]
			if( !is.function(f <- args[[ifun]]) )
				stop("Invalid scoring method definition: first unnamed argument must be a function [", class(f), ']')
			args <- c(list(f), args[-ifun])
			setPackageRegistryEntry('score', key=method, object=object, algorithm=args)
			f
		} else GEDscore_registry$delete_entry(method)
	}
}


#' Marker Scoring Method: Abbas et al. (2009)
#' 
#' \code{markerScoreAbbas} implements the scoring/selection method proposed 
#' by \cite{Abbas2009}, to select marker genes from pure cell type samples.
#' 
#' The method \sQuote{Abbas} uses a t-test approach, so that the data is 
#' assumed to contain at least 2 pure samples per cell-type. 
#' It implements the method from \cite{Abbas2009}:
#' 
#' "[...] top differentially expressed (based on 95% fold change confidence intervals 
#' from Student's T-test) probesets were determined by comparing each
#' probe's highest-expressed group with the next highest-expressed
#' group in order to find probesets that are good markers for each cell
#' population. This step was repeated with comparison between the
#' top group and the third-highest group in order to also include
#' probesets that were strong markers for two cell populations."
#' 
#' For each gene, the highest-expressing cell type is determined by ordering them 
#' by mean expression.
#' Comparisons and p-value computations are performed using the fast t-test 
#' implementation from \code{\link{rowttests}} in the \pkg{genefilter} package.
#' 
#' @inheritParams extractMarkers
#' @param statistic statistic to use as a score.
#' The method computes the following quantities for each comparison 
#' between the most expressing cell type (highest within-group mean) with the 
#' i-th most expressing other cell type:
#' \describe{
#' \item{p.value<i>}{p-value (t-test).}
#' \item{dm<i>}{difference in means.}
#' \item{statistic<i>}{statistic value (t statistic).}
#' \item{dmM<i>}{minimum difference = min(top) - max(top<i>).}
#' \item{fold<i>}{fold change in means.}
#' \item{mMfold<i>}{minimum fold change = min(top)/max(top<i>).}
#' \item{top<i>}{group index of the i-th most expressing cell type.
#' It is relative to the levels of \code{data}.}
#' }
#' @param ntop the number of groups for which statistics should be 
#' computed (>= 2).
#' If \code{ntop=2}, only statistics between the highest and 
#' second-highest expressing cell-type are computed.
#' If \code{ntop=3}, statistics between the highest and third-highest are also 
#' computed, and so on.
#' Use \code{ntop=Inf} to compare the highest-expressing group to all other groups.
#'  
#' @param log a logical specifying if the should the data be log-transformed before computing the 
#' p-values.
#' Default is to log-transform the data if is not already in log-scale, which 
#' is determined by the function \code{\link{is_logscale}}.
#' 
#' @param lbase log base to used if the data is log-transformed. 
#' @param vsall logical used when \code{ntop=2} that indicates if the 
#' comparison should be carried out between the highest-expressing cell-type
#' and the rest of all other cell-types.
#'  
#' @return \item{extractMarkers}{a numeric matrix with the following named columns:
#' \itemize{
#' \item{top}{the index of the highest-expressing cell-type, as defined by the 
#' levels of the \code{factor} derived from argument \code{group}.}
#' \item{p.value2}{ p-value from the comparison between the highest-expressing 
#' (top) and second highest-expressing (second top) cell-type.}
#' \item{dm2}{difference in mean expression between the top and second top cell-type}
#' \item{dmM2}{minimum difference of expression between the top and second top cell-type}
#' \item{fold2}{fold change mean expression between the top and second top cell-type}
#' \item{mMfold2}{minimum fold change expression between the top and second top cell-type}
#' \item{dm3, mdM3, fold3, mMfold3}{ same quantitities computes between the 
#' top and third-top cell-type.}
#' }
#' 
#' The result matrix has an attribute \code{'types'} that contains the levels 
#' of the original (or converted) factor \code{group}.}
#'  
#' @family markerScore
#' @importFrom genefilter rowttests
#' @export
markerScoreAbbas <- markerScoreMethod('Abbas', object='matrix', 
	function(object, data, statistic='p.value', ntop=2, log=!is_logscale(object), lbase=2, vsall=FALSE){
				
		if( missing(data) || is.null(data) )
			stop("Missing argument `group` for marker scoring method 'Abbas': cell type groups are required")
			
		group <- data
		nl <- nlevels(group)
		if( nl < 2 )
			stop("extractMarkers - Invalid argument: `data` must have more than 1 level.")
		if( !missing(ntop) && is.finite(ntop) && ntop > nl ){
			warning("NOTE: Number of comparisons [", ntop, "] exceeds number of groups [", nl, "].")
		}
		ntop <- max(min(ntop, nl), 2)
		
		# log-transform
		if( is.null(log) ){
			log <- !is_logscale(object)
		}else if( isNumber(log) ){
			lbase <- log
			log <- TRUE
		}
		if( isTRUE(log) )
			object <- log_transform(object, lbase)
		#
		
		# compute t-test between top-group and second top-group
		group_idx <- 1:nlevels(group)
		group_list <- split(seq_along(group), group) 
		
#  		objectVal <- rowMaxsBy(object, group)
# 		objectVal <- t(sapply(1:nrow(object), function(i){
# 			sapply(group_list, function(a) min(object[i, a])/max(object[i, -a]))
# 		}))
# 		print(sum(objectVal>1))
		objectVal <- rowMeansBy(object, group)
		
		wm <- t(apply(objectVal, 1L, 
			function(x){
				# first top
				ig <- which.max(x)
				# second top
				res <- c(ig, group_idx[-ig][which.max(x[-ig])])
				# other top groups
				if( ntop > 2 ){
					for( i in 3:ntop ){
						res <- c(res, group_idx[-res][which.max(x[-res])])
					}
				}
				res
			}
		))
		
		# local function to compute t-tests between top group and another top group
		.compute_t <- function(j, i=1){
			# build list of all combinations
			l <- paste(wm[,i], wm[,j])
			l <- factor(l, levels=unique(l))
			top <- lapply(levels(l), function(x){
						i <- as.integer(strsplit(x, ' ', fixed=TRUE)[[1]])
						list(group_list[[i[1L]]], group_list[[i[2L]]]) 
					})
			idx <- split(1:nrow(object), l)
			var.equal <- TRUE
			tt <- mapply(function(i, g){
						#			res <- fastT(object[i, c(g[[1L]], g[[2L]]), drop=FALSE]
						#				, seq_along(g[[1L]]) 
						#				, length(g[[1L]]) + seq_along(g[[2L]]), var.equal=var.equal)
            if( ntop == 2 && vsall ) g[[2L]] <- (1:ncol(object))[-g[[1L]]]
						res <- as.matrix(rowttests(object[i, c(g[[1L]], g[[2L]]), drop=FALSE]
										, factor(c(rep(1, length(g[[1L]])), rep(2, length(g[[2L]]))))))
						res <- res[, c(3,2,1), drop=FALSE]
						# add other stats
						x1 <- object[i, g[[1L]], drop=FALSE]
						x2 <- object[i, g[[2L]], drop=FALSE]
						m1 <- rowMeans(x1); m2 <- rowMeans(x2)
						min1 <- rowMins(x1);
						max2 <- rowMaxs(x2)
						
						cbind(res, dmM= min1 - max2, fold=m1/m2, mMfold= min1/max2)
					}, idx, top, SIMPLIFY=FALSE)
			
			tt <- do.call('rbind', tt)
			# order as in the original matrix 
			tt <- tt[order(unlist(idx)), ]
			# append top index to colnames
			cn <- if( i==1 ) j else str_c(i, '_', j)
			colnames(tt) <- paste(colnames(tt), cn, sep='')
			tt
		}
		
		# compute stats against 2nd top group
		res <- cbind(top=wm[, 1], .compute_t(2), wm[, 2])
		colnames(res)[ncol(res)] <- 'top2'
		# compute stats against 3rd top group
		if( ntop > 2 ){
			for( j in 3:ntop ){
				res <- cbind(res, .compute_t(j), wm[, j])
				colnames(res)[ncol(res)] <- str_c('top', j)
			}
		}
		# compute min p-value
		res <- cbind(res, p.value=rowMins(res[, grep("^p\\.value", colnames(res)), drop=FALSE]))
		
		# check for NA values
	    if( anyMissing(pv <- res[, 'p.value']) ){
	      rowid <- which(is.na(pv))
	      if( !is.null(rn <- rownames(res)) ) rowid <- rn[rowid]
	      warning("Some of the p-values could not be computed: NA values introduced for ", length(rowid), " feature(s) [", str_out(rowid), "].")
	    }
		
	    # copy selected score statistic column in second
		res <- cbind(res[,1L,drop=FALSE], res[, statistic], res[,-1L,drop=FALSE])
		colnames(res)[2] <- statistic
		
		# result
		return(res)
	}
)

#' Marker scoring method: Max Expression 
#' 
#' This method scores genes according with their highest expression in a set of pure
#' samples.
#' 
#' Genes are assigned to their respective maximum expressing column/value and returned
#' as a \code{MarkerList} object.    
#' 
#' @inheritParams extractMarkers
#' 
#' @family markerScore
#' @export
#' @examples
#' 
#' # random data matrix with rownames (needed)
#' x <- rmatrix(20, 3)
#' rownames(x) <- paste("g", 1:nrow(x), sep='_')
#' 
#' # compute score
#' markerScoreMaxcol(x)
#' # same thing with auto-ordering through the main interface function 
#' extractMarkers(x, method='maxcol')
#' 
markerScoreMaxcol <- markerScoreMethod('maxcol', object='matrix', decreasing=TRUE, 
	function(object){

		# get macimum contributing column
		imax <- max.col(object, ties.method='first')
		if( is.null(colnames(object)) ) colnames(object) <- 1:ncol(object)
		# assign each row to its col max value
		m <- unlist(mapply(function(i, j) object[i, j], seq(nrow(object)), imax))
		# require rownames
		if( is.null(nm <- rownames(object)) )
			stop("Cannot create marker list: input matrix data must have rownames to associate values to rows.")
		m <- setNames(m, nm)
		# return as a MarkerList object
		types <- factor(colnames(object)[imax], levels=colnames(object))
		MarkerList(m, names=types)
	}
)


#' Marker scoring method: SCOREM 
#' 
#' This method applies an approach called \emph{statistical consolidation of redundant 
#' expression measures} (SCOREM) from \cite{Schneider2011} whitin each cell type, 
#' in order to select groups of markers that have consistent expression in 
#' mixed samples.
#' 
#' The SCOREM approach allows sets of marker genes to be optimised with respect to 
#' the dataset one intends to use them with.
#' It uses a sub-graph detection algorithm, on a graph defined from
#' the rank correlation matrix of marker expression profiles.
#' 
#' \strong{IMPORTANT: this feature is still experimental and in development.}
#' 
#' @inheritParams extractMarkers
#' @param minsd minimum standard deviation of expression marker profiles
#' @param alpha numerical value between 0 and 1, indicating the level of
#' significance required for correlations to be considered significant.
#' It is a p-value threshold: the lower the value the more markers  
#' must be correlated to be called in the same sub-group.
#' @param ... extra arguments not used in \code{markerScoreScorem}, 
#' but passed to \code{\link{selectMarkers.MarkerList}} in \code{selectMarkers}.
#' @param verbose logical that toggles log messages (passed to \code{\link{subset,MarkerList-method}}), 
#' which can be useful to track any ID conversion.
#' 
#' @family markerScore
#' @export
#' @examples
#' 
#' # random data matrix and marker list
#' x <- rmix(3, 100, 20)
#' m <- getMarkers(x)
#' 
#' # compute SCOREM groups: should keep all markers together
#' # or remove weak markers 
#' s <- markerScoreScorem(m, x)
#' str(s)
#' # add some non marker genes
#' m2 <- m
#' m2[[1]] <- as.integer(c(m2[[1]], c(40, 50, 60)))
#' s <- markerScoreScorem(m2, x)
#' str(s[[1]])
#' 
#' # with group selection and reordering through the main interface function 
#' extractMarkers(m, x, method='scorem')
#' extractMarkers(m2, x, method='scorem')
#' 
#' @demo Data-driven marker selection (experimental)
#' # Using the SCOREM approach one can filter marker genes that have consisten 
#' # whose expression values in a given datasets.
#' # This is because the notion of marker can be very dependent on the data, biological 
#' # conditions or platform used to define them.
#' 
#' # load kidney transplant dataset from Shen-Orr et al. (2010)
#' eset <- ExpressionMix('GSE20300', verbose=TRUE)
#' 
#' # load HaemAtlas marker gene list (on Illumina)
#' ml <- MarkerList('HaemAtlas')
#' ml
#' # convert marker IDs to Affy IDs from dataset
#' # using stringent one to one mapping
#' mla <- convertIDs(ml, eset, method='1:1', verbose=TRUE)
#' summary(mla)
#' 
#' # plot expression profile of each set of markers (only the first 10)
#' profplot(mla[, 1:10], eset, split=TRUE, lab='')
#' 
#' # filter out using SCOREM
#' sml <- extractMarkers(mla, eset, method='scorem', alpha=.005)
#' # plot selected markers
#' profplot(sml, eset, split=TRUE, lab='')
#' 
#' 
markerScoreScorem <- markerScoreMethod('scorem', object='MarkerList', decreasing=TRUE, 
	function(object, data, minsd=0.1, alpha=.01, ..., verbose=FALSE){
		
		# put MarkerList in `object`
		if( isMarkerList(data) ){
			swap <- object
			object <- data
			data <- swap
		}
		if( !isMarkerList(object) )
			stop("Invalid marker list object [", class(object), ']')
		
		# subset marker list to match features in expression data (potential ID conversion here)
		object <- subset(object, data, verbose=verbose)
		# make scorem groups
		res <- scoremGroups(object, data, minsd=minsd, alpha=alpha, verbose=verbose)
    # pass on new annotation as an attribute
    attr(res, 'annotation') <- annotation(object)
    # return result
    res
	}
)

#' The method \code{selectMarkers} selects, within each cell type separately,
#' the group of markers that have the more consistent expression values, 
#' according to the criterion specified in argument \code{select}:
#' \describe{
#' \item{\code{'score'}}{selects group(s) with the highest SCOREM score (i.e. mean spearman correlation).}
#' \item{\code{'size'}}{selects group(s) of maximum size.}
#' }
#' 
#' If multiple groups pass the criterion, then the group with the lowest
#' median gene expression gene-wise standard deviation is selected.
#' This is to select the group of markers whose expression values are 
#' the more similar.
#' 
#' Before returning the \code{MarkerList} object, each selected marker is assigned a 
#' score defined as its minimum correlation with other markers in its group.
#' 
#' @inheritParams selectMarkers
#' @param mscore correlation method used to compute a score for each marker of the selected
#' groups.
#' 
#' @S3method selectMarkers markerScore_scorem
#' @rdname markerScoreScorem
selectMarkers.markerScore_scorem <- function(x, data, mscore='pearson'
											, statistic=c('size', 'score'), ..., .object=NULL){
	
	# match selection criterion specification
	statistic <- match.arg(statistic)
	
	# select group with the best correlation
	res <- sapply(x, function(g){
				
		# get sub-group(s) with max score
		lg <- sapply(g, function(x) length(x$ids))		
		# a priori look for non-trivial sub-groups
		subs <- g[lg > 1]
		if( length(subs) > 0L ) g <- subs
		else if( length(g) > 1L ){
			# choose first one
			# TODO: improve this, e.g., choose most stair-like in pure samples
			g <- g[1L]
		}
		s <- sapply(g, '[[', 'score')
		
		# compute criterion
		if( statistic == 'score' ){ # max SCOREM score
			mg <- g[s == max(s)]
		}else if( statistic == 'size'){ # max size
			lg <- sapply(g, function(x) length(x$ids))	
			mg <- g[lg == max(lg)]
		}
		
		# take min variance group if multiple hits
		if( length(mg) > 1L ){ 
			mvar <- sapply(mg, function(x){
				median(colSds(exprs(data)[x$ids, , drop=FALSE]))
			})
			mg <- mg[which.min(mvar)]			
		}else if( length(mg) == 0L ){ # only singletons
		  warning("no marker sub-group meet criterion")
		  g[s == max(s)]
		}
		mg <- mg[[1L]]
		
		# extract ids
		ids <- mg$ids
		m <- exprs(data)[ids,, drop=FALSE]
		
		# assign intra-group correlation score to each marker
		if( length(ids) > 1 ){
			co <- cor(t(m), method=mscore)
			setNames(rowMins(co, na.rm=TRUE), ids)
		} else if( length(ids) == 1L ) setNames(1, ids)
		else numeric()
	}, simplify=FALSE)

	if( is.null(.object) ) .object <- MarkerList(res) 
	else geneIds(.object) <- res
  
	# set new annotation if necessary
  if( hasAnnotation(x) && !identical(getAnnotation(x), getAnnotation(.object)) ){
    annotation(.object) <- NULL
    annotation(.object) <- getAnnotation(x)
  }
	
	# call generic selection method
	selectMarkers(.object, ...)
}

#' Marker Scoring Method: Tukey Honest Significant Difference
#' 
#' This scoring method \code{markerScoreHSD} performs pairwise 
#' comparison between groups of pure samples, and scores each 
#' comparison using Tukey's Honest Significant Difference 
#' p-values (see \code{\link{TukeyHSD}}).
#' 
#' The scores are returned in a matrix, with features in rows 
#' and cell types in column, which contains the HSD p-values 
#' corresponding to the comparisons between the most expressing 
#' cell type and  other cell types.
#' Each row contains an NA value that identifies the column 
#' corresponding to the associated feature's most expressing cell type. 
#' 
#' Features whose expression is not consistently higher in one cell type 
#' than in any other cell type are discarded.
#' 
#' @seealso \code{\link{TukeyHSD}}
#' @inheritParams markerScoreAbbas
#' @param verbose verbosity level, usually \code{TRUE} or \code{FALSE}.
#' @family markerScore
#' @export
#' @examples 
#' 
#' # generate data from pure cell type samples
#' x <- rpure(3)
#' x
#' aheatmap(x, annCol=TRUE)
#' 
#' # extract markers
#' ml <- extractMarkers(x, x$CellType, method='HSD')
#' # check score/p-value distribution
#' hist(ml)
#' # plot most significant ones
#' profplot(ml[ml < 0.0001], x, split=TRUE)
#' 
markerScoreHSD <- markerScoreMethod('HSD', object='matrix', decreasing=FALSE,  
	function(object, data, log=!is_logscale(object), lbase=2, verbose=FALSE){
		
		if( missing(data) || is.null(data) )
			stop("Missing argument `group` for marker scoring method 'HSD': cell type groups are required")
		
		nl <- nlevels(data)
		if( nl < 2 )
			stop("extractMarkers - Invalid argument: `data` must have more than 1 level.")
		
		# log-transform
		if( is.null(log) ){
			log <- !is_logscale(object)
		}else if( isNumber(log) ){
			lbase <- log
			log <- TRUE
		}
		if( isTRUE(log) )
			object <- log_transform(object, lbase)
		#
		
		# pre-compute indexes and sign correction to extract Tukey's test results
		# for each cell type cell-type
		ctnames <- levels(data)
		fit <- aov(object[1,] ~ data)
		t <- TukeyHSD(fit)
		.pair <- function(ct){
			others <- ctnames[ctnames != ct]
			c(paste(ct, others, sep='-'), paste(others, ct, sep='-'))
		}
		ix <- sapply(ctnames, function(ct) rownames(t$data) %in% .pair(ct))
		# check result
		cs <- colSums(ix)
		if( !all(cs ==  cs[1]) )
			stop("Unexpected error: cell-type pairwise comparison index names not as expected")
		sgn <- sapply(ctnames, function(ct){
			others <- ctnames[ctnames != ct]
			last <- rownames(t$data) %in% paste(others, ct, sep='-')
			ifelse(last, -1, 1)
		})
		
		# compute difference from each cell type: must all be positive
		# init result matrix
		pvalues <- NA_Matrix(nrow(object), length(ctnames))
		colnames(pvalues) <- ctnames
		# do compute
		idx <- 0
		if( verbose ){
			message("Computing HSD statistics")
			pgr <- iterCount(nrow(object), title='Features')
		}
		mn <- lapply(1:nrow(object), function(f){
			if( verbose ) pgr(f)
			x <- object[f,]
			# fit anova
			fit <- aov(x ~ data)
			# compute Tukey's test
			t <- TukeyHSD(fit)
			hsd <- t$data
			# reverse differences where needed and flag inconsistently expressed features
			sapply(1:ncol(ix), function(k){
				i <- ix[, k]
				t.ct <- hsd[i, , drop=FALSE]
				# only consider features that are constistently over expressed
				diff <- sgn[i, k] * t.ct[, 1L]
				if( any(diff <= 0 | t.ct[, 4L] < 0) ) NA
				else{
					idx <<- idx+1
					pvalues[idx, -k] <<- t.ct[, 4L]
					f
				}
			})
		})
	
		# drop unecessary data
		pvalues <- head(pvalues, idx)
		
		mn <- unlist(mn)
		mn <- mn[!is.na(mn)]
		stopifnot( length(mn) == nrow(pvalues) )
		# check unexpected multiple cell-types
		if( anyDuplicated(mn) ){
			warning("Skipped markers ", str_out(mn[duplicated(mn)]), ": specific to multiple cell types ["
							, "]")
			nodup <- !mn %in% mn[duplicated(mn)]
			pvalues <- pvalues[nodup, , drop=FALSE]
			mn <- mn[nodup]			
		}
		# use rownames if any
		if( !is.null(rownames(object)) ) rownames(pvalues) <- rownames(object)[mn]
		else rownames(pvalues) <- mn
		
		# return p-values
		pvalues
	}
)

#' The method \code{selectMarkersmarkerScore_HSD} selects, 
#' within each cell type separately, the markers with the lower 
#' aggregated p-value for Tukey's HSD.
#' The default aggregation method is to compute the maximum HSD p-value.
#' 
#' @inheritParams selectMarkers
#' @param statistic method used to aggregate p-values of the pairwise comparisons 
#' with each of the other cell types into a single numeric score.
#' 
#' @S3method selectMarkers markerScore_HSD
#' @rdname markerScoreHSD
selectMarkers.markerScore_HSD <- function(x, data, statistic=max, ...){

	# compute score
	s <- apply(x, 1L, function(x){
		statistic(x[!is.na(x)])
	})
	# split cell types
	cti <- apply(x, 1L, function(x) which(is.na(x)))
	ct <- colnames(x)[cti]
	ml <- MarkerList(split(s, ct))
	# order/select based on threshold
	selectMarkers(ml, ...)
}

