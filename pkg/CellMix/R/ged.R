# Main interface methods/functions for gene expression deconvolution 
# 
# Author: Renaud Gaujoux
# Creation: 19 Jan 2012
###############################################################################

#' @include AllClasses.R
#' @include GEDStrategy-class.R
NULL

#' Main Interface for Gene Expression Deconvolution Methods
#'
#' @param object global gene expression matrix-like data object 
#' (e.g., \code{matrix} or \code{ExpressionSet})
#' @param x input data used by the algorithm to deconvolve global gene expression.
#' @param method specification of a gene expression deconvolution method.
#' If missing, a default suitable algorithm is selected, based on the type and/or 
#' dimensions of \var{x} and \var{data}.
#' @param ... extra arguments to allow extension, most of which are passed down 
#' to the deconvolution algorithm itself. 
#' 
#' @inline
#' @export
#' @examples
#' 
#' # random global expression data: 3 cell types, 20 samples, 100 features
#' X <- rmix(3, 100, 20, markers=5)
#' dim(X)
#' # extract signature/proportion/markers
#' sig <- basis(X)
#' prop <- coef(X)
#' m <- getMarkers(X)
#' summary(m)
#' 
setGeneric('ged', function(object, x, method, ...) standardGeneric('ged') )
#' Default \code{ged} method apply the auto-selection scheme for determining 
#' which method is suitable for the type of input data.
#'
#' @param data optional data, typically a marker list, specified in a format 
#' that is supported by the factory function \code{\link{MarkerList}}.
#' @param maxIter maximum number of iterations to perform.
#' If \code{method} is missing, the value of this argument can influence which 
#' method is selected. 
#' See section \emph{Details}.
#' @param verbose logical that toggles verbosity.
#' A number (integer) can be passed to specify the verbosity level (the higher the more messages 
#' are output).
#' Passing \code{verbose=Inf} toggles debug mode (all messages).
#' Note that because it appears after \code{...} it must be fully named.
#' @param wrap logical that indicates the result returned by the method should be wrapped 
#' into an NMF object.
#' If \code{FALSE}, it is returned as is.
#' Note that because it appears after \code{...} it must be fully named.
#' @param dryrun logical that indicates if one should make a fake run, i.e., without fitting any data,
#' but showing as much details about the run as possible (used or inferred parameters, etc..).  
#' Note that because it appears after \code{...} it must be fully named.
#'  
setMethod('ged', signature(object='MatrixData', method='GEDStrategy')
, function(object, x, method, data=NULL, maxIter=1L, ...
			, verbose=cellmix.getOption('verbose'), wrap=TRUE, dryrun=FALSE){
	
	gedopt <- .cellmix.options(wrap=wrap)
	on.exit({ .cellmix.options(gedopt) }, add=TRUE)
	
	# set verbosity level
	if( !missing(verbose) || verbose ){
		ol <- lverbose(verbose)
		opt.nmf <- nmf.options(verbose=max(verbose-1, 0))
		on.exit({ lverbose(ol); nmf.options(opt.nmf) }, add=TRUE)
	}
	
	# load markers as a MarkerList if present
#	if( !missing(data) && is.list(data) && !is.data.frame(data) && !is.matrix(data) ){
#		lmessage(2L, "ged - Loading MarkerList object from argument `data` [<", class(data), ">]")
#		data <- MarkerList(data)
#	}
	
	if( missing(x) || is.null(x) ){
		if( !isMarkerList(data) )
			stop("ged - Invalid parameters: cannot infer number of cell types with both `x` and `data` both missing or NULL.")
		x <- length(data)
	}
	
	# extract method name
	methodName <- name(method)
	
	# force using gedProportion pipeline if necessary
	# send to gedProportions if suitable and not already within a call to it
	if( !.in_gedProportion() && isMatrixData(x) && isRequired('Basis', method) ){
		if( verbose ) logger.options(autoindent=FALSE)
		return( gedProportions(object, x, method = method, ..., data = data, maxIter = maxIter
								, verbose = verbose, wrap = wrap, dryrun = dryrun) )
	}
	#
	
	# inform the user of the selected algorithm
	vmessage("Using ged algorithm: ", dQuote(name(method)))
	
	# details on the strategy
	if( verbose > 3L ){
		vmessage("## GED strategy:")
		s <- capture.output(show(method))
		vmessage(paste('#', s, collapse="\n"))
		vmessage('##')
	}
	
	# use default maxIter value if necessary
	.call <- as.list(as.call(list(quote(method@algorithm), quote(object), quote(x))))
	if( !missing(maxIter) && hasFormals(method@algorithm, 'maxIter') )
		.call <- c(.call, maxIter=quote(maxIter))

	# limit markers if necessary
	if( isMarkerList(data) && isNumber(x) ){
		rank <- x
		# adjust length of markers if necessary
		if( rank > length(data) ){
			warning("ged - Number of cell types greater than number of marker types: some components won't use markers")
			# extend with empty cell markers
			n.extra <- rank-length(data)
			for(i in 1:n.extra ){
				data[[str_c("Unknown_", length(data)+i)]] <- character()	
			}
		}else if( rank < length(data) ){
			lmessage(2L, "ged - Using only first ", rank, " out of ", length(data), " marker types out of ")
			data <- head(data, rank)
		}
	}
	# pass data to the strategy's algorithm
	if( !is.null(data) ){
		.call <- c(.call, data=quote(data))
	}
	# add default arguments
	if( length(method@defaults) )
		.call <- expand_list(.call, method@defaults, .exact=FALSE)
	.call <- c(.call, quote(...))
	.call <- as.call(.call)
#	print(.call)
	# init result before fitting so that the intial random seed state is correctly stored
	nmfres <- NMF:::NMFfit(method=methodName)
	
	# DO FIT
#	print(.call)
	t <- system.time(res <- eval(.call))
	##

	vmessage("Timing:")
	if( verbose ) show(t)
	
	# cleanup memory
	gc()
	# return raw result if requested
	if( !wrap ) return(res)
	
	vmessage("GED final wrap up ... ", appendLF=FALSE)
	# wrap into an NMF fit object
	if( !isNMFfit(res) ){
	
		# extract expression data if necessary
		if( isExpressionSet(res) ) res <- exprs(res)
		# wrap 
		if( is.matrix(res) ){
			if( ncol(res) == ncol(object) ){# proportions
				# force dimnames
				colnames(res) <- sampleNames(object)
				res <- nmfModel(H=res)
			}else if( nrow(res) == nrow(object) ){# signatures
				# force dimnames
				rownames(res) <- featureNames(object)
				res <- nmfModel(W=res)
			}
		}
		# set model in fit object
		if( is.nmf(res) ){
			fit(nmfres) <- res
			nmfres@runtime <- t
			res <- nmfres
		}
	}else{
		# enforce ged algorithm name to overwrite possible NMF algorithm name
		algorithm(res) <- methodName
	}
	
	# force dimnames
	if( !is.list(res) ){
		if( nrow(res) && is.null(rownames(res)) && nrow(res) == nrow(object)){
			rownames(res) <- featureNames(object)
		}
		if( ncol(res) && is.null(colnames(res)) && ncol(res) == ncol(object) ){
			colnames(res) <- sampleNames(object)
		}
	}else{
	}
	#

	# add S3 class to raw result
	# set class on element `rawfit`
	s3cl <- str_c('gedfit_', methodName)
	if( !is.null(res$rawfit) ){
		class(res$rawfit) <- s3cl 
	}else if( is.list(res) && !isS4(res) ){
		class(res) <- c(s3cl, class(res))
	}
	#

	# compute statistics with csSAM
#	csfit <- ged(object, coef(res), 'csSAM')
#	if( !hasBasis(res) )
#		res <- csfit
#	#
	
	vmessage('OK')
	# return result
	res
}
)

#' Applies a deconvolution algorithm registered in the \pkg{CellMix} 
#' registry.
setMethod('ged', signature(object='MatrixData', method='ANY'),
	function(object, x, method, ...){
		# create strategy based on input
		method <- gedAlgorithm(method)
		ged(object=object, x=x, method=method, ...)
	}
)


#' This method deconvolves the target expression matrix using 
#' the algorithm selected by the automatic selection strategy 
#' implemented by \code{\link{selectGEDMethod}}, which choose 
#' a suitable algorithm whose data requirements match the 
#' provided input data (i.e. arguments \code{x} and optionally \code{data}).
#' See some more details in the \code{\link{selectGEDMethod}}
#' man page.
#' 
#' @seealso \code{\link{selectGEDMethod}}
#'  
#' @examples
#' #--------------------------------------------
#' # Automatic selection of a suitable algorithm
#' #--------------------------------------------
#' # expression data only: fastdeconf
#' res <- ged(X, 3)
#' # with markers only, non iterative: qprog
#' res <- ged(X, m)
#' # with markers only, iterative: ssKL
#' res <- ged(X, m, maxIter=5)
#' # with signatures: lsfit
#' res <- ged(X, sig)
#' # with proportions: csSAM
#' res <- ged(X, prop)
#' # with proportions, iterative: DSection
#' if( require.quiet(RcppOctave) ){
#' 	res <- ged(X, prop, maxIter=5)
#' }
#' 
setMethod('ged', signature(object='MatrixData', method='missing'),
	function(object, x, method, ...){
		ged_call <- selectGEDMethod(object, x=x, ..., call. = TRUE, quiet.=FALSE)
		eval(ged_call)
	}
)

#' Automatic Selection of Gene Expression Deconvolution Algorithms
#' 
#' Implements a simple automatic selection strategy that chooses
#' a suitable deconvolution method given the input data.
#' 
#' The selection aims at finding an algorithm that is able to 
#' perform deconvolution based on the provided input data.
#' The strategy is to choose amongst the possible algorithms 
#' available from the \pkg{CellMix} \emph{built-in} registry, 
#' according to their respective data requirements.
#' 
#' Essentially the choice of algorithms made based on the dimensions
#' of the target expression data \code{object} and the dimensions or 
#' type of the input data in \code{x} and \code{data}. 
#' 
#' Currently, the pipeline does not attempt is made to 
#' choose the "best" algorithm, which would be the one 
#' that would return the most accurate results (proportions 
#' or cell-specific signatures/differences) for the given data
#' setting.
#' 
#' The selected algorithm is indeed chosen as to be  
#' \emph{applicable} to the input data.
#' When possible, however, a state of the art algorithm or the 
#' most currently used algorithm is selected.
#' 
#' @inheritParams ged
#' @param call. logical that indicates if one should return 
#' the suitable call to \code{\link{ged}} (\code{TRUE}), or
#' just the name of the selected method.
#' @param quiet. logical that indicates if an error should be 
#' thrown if no algorithm able to fit the input data is found, 
#' or simply return \code{NULL}.
#' 
#' If explicitly set to \code{FALSE}, then a note is displayed,  
#' showing the selected algorithm and a quick justification for  
#' the choice.
#' 
#' @export 
#' @examples
#' # ged methods requirements
#' selectGEDMethod()
#' # generate mixed expression data
#' x <- rmix(3, 100, 20)
#' dim(x)
#' sig <- basis(x)
#' prop <- coef(x)
#' ml <- getMarkers(x)
#' 
#' # one need at least the number of cell types 
#' try( selectGEDMethod(x) )
#' selectGEDMethod(x, 3)
#' # from signature basis matrix
#' selectGEDMethod(x, sig)
#' selectGEDMethod(x, sig, quiet.=FALSE)
#' # from cell proportion matrix
#' selectGEDMethod(x, prop)
#' # from cell proportion matrix with multiple iterations
#' selectGEDMethod(x, prop, maxIter=10, quiet.=FALSE)
#' # from cell proportion matrix with markers
#' selectGEDMethod(x, prop, data=ml, quiet.=FALSE)
#' # from marker genes
#' selectGEDMethod(x, ml)
#' # from marker genes with multiple iterations
#' selectGEDMethod(x, ml, maxIter=10, quiet.=FALSE)
#' 
selectGEDMethod <- function(object, x=NULL, data=NULL, maxIter=1L, ..., call. = FALSE, quiet. = FALSE){
		
	# return ged method requirements
	if( nargs() == 0L ) return( gedRequired() )
	
	# check parameters
	if( !isMatrixData(object) ){
		stop('Target data must be a matrix-like object [', class(object), ']')
	}
	# no note if quiet mode is required
	if( missing(quiet.) ) shownote <- FALSE
	else shownote <- !quiet.
	# build common error message header
	errmsg <- paste0("Could not find a suitable GED method for input: "
			, str_out(c(object=class(object)
						, x=if( missing(x) ) 'missing' else class(x)
						, data=if( missing(data) ) 'missing' else class(data), maxIter=maxIter), Inf, use.names=TRUE), ".")
	algos <- paste0("\n  Try explicitly specifying a GED algorithm: ", str_out(gedAlgorithm(), Inf))
	
	rank <- x
	# extract dimensions
	dt <- dim(object) # target
	dx <- dim(x) # input
	
	reason <- NULL
	if( is.null(dx) ){ # not a matrix-like object
		
		if( isMarkerList(x) ){ # select a method that uses marker genes
			if( maxIter <= 1L ){
				reason <- "markers + maxIter=1"
				method <- 'DSA'
			} else{
				reason <- "markers + maxIter>1"
				method <- 'ssKL'
			}
		}else if( is.numeric(x) ){ # x a numeric vector			
			rank <- x
			if( is.null(data) ){ # use deconf
				reason <- "no prior data"
				method <- 'deconf'
			}
			else if( isMarkerList(data) ){
				if( maxIter <= 1L ){
					reason <- "markers + maxIter=1"
					method <- 'DSA'
				}
				else{
					reason <- "markers + maxIter>1"
					method <- 'ssKL'
				}
			}else{ # not supported
				if( quiet. ) return() # return NULL
				stop(errmsg, "\n  Supported selection formats for 'data' when 'x' is not a matrix-like object are: numeric vector or a MarkerList object."
							, algos)
			}
		}else{ # not supported
			if( quiet. ) return() # return NULL
			stop(errmsg, "\n  Supported selection formats for 'x' are: missing, matrix-like objects, numeric vectors or MarkerList objects."
					, algos)
		}
	}else if( (dx[1] == dt[1]) || (dx[1] > dx[2]) ){ # x is a basis matrix: use lsfit
		# - same number of rows in x and target
		# - more rows than columns in x
		rank <- dx[2]
		reason <- "signature data"
		method <- 'lsfit'
	}else if( (dx[2] == dt[2]) || (dx[2] > dx[1]) ){# x is a proportion matrix
		rank <- dx[1]
		reason <- "proportion data"
		if( maxIter <= 1L ){ # Non-iterative algorithm
			reason <- str_c(reason, ' + maxIter=1')
			
			if( is.null(data) || is.vector(data) || is.factor(data) ){
				if( is.vector(data) || is.factor(data) )
					reason <- str_c(reason, ' + groups')
				method <- 'csSAM'

			}else if( isMarkerList(data) ){
				reason <- str_c(reason, ' + markers')
				method <- 'cs-qprog'

			}else{
				if( quiet. ) return() # return NULL
				stop(errmsg, "\n  Detected 'x' as a proportion matrix and single iteration fit.\n"
							, "  Supported selection formats for 'data' in this case are : NULL, vector, factor or MarkerList objects."
							, algos)
				
			}
		}else{ # Iterative algorithm
			reason <- str_c(reason, ' + maxIter>1')
			
			if( is.null(data) || is.vector(data) || is.factor(data) ){
				reason <- str_c(reason, ' + groups')
				method <- 'DSection'
			}else{
				if( quiet. ) return() # return NULL
				stop(errmsg, "\n  Detected 'x' as a proportion matrix and multiple iteration fit.\n"
						, "  Supported selection formats for 'data' in this case are : NULL, vector or factor objects."
						, algos)
			}
		}
	}else{
		if( quiet. ) return() # return NULL
		stop(errmsg, "\n  Detected 'x' as matrix-like data.\n"
				, "\n  Supported selection formats for 'data' in this case are : NULL, vector or factor objects."
				, algos)
	}
	
	# inform the user of the selected algorithm
	vmessage(if( shownote ) "# NOTE - "
			, "Selected GED algorithm: ", dQuote(method), ' [', reason, ']', force=shownote)
	
	# return selected method or call if requested
	if( call. ){
		cl <- match.call()
		cl['call.'] <- NULL
		cl['quiet.'] <- NULL
		cl['method'] <- method
		cl[[1L]] <- as.name('ged')
		as.call(cl)
	}else method	
}

## #' Applies an gene expression algorithm from a list of markers only.
## #' 
## #' When \var{data} is not specified, it is a shortcut for 
## #' \code{ged(object, method=method, data=x, ...}.
## #'  
#setMethod('ged', signature(object='ANY', x='MarkerList'),
#	function(object, x, method, data=NULL, ...){
#		
#		if( missing(method) ) method <- NULL
#		# use markers passed in x if no other data is otherwise specified
#		if( is.null(data) ){
#			ged(object, method=method, data=x, ...)
#		}else{
#			callNextMethod()
#		}
#
#	}
#)

#' Applies a custom gene expression algorithm.
#' 
#' @param name optional name for the custom algorithm.
setMethod('ged', signature(method='function'),
	function(object, x, method, ..., name=NULL){
		
		method <- GEDStrategy(method, name=name)
		if( missing(x) ) x <- NULL
		callGeneric(object=object, x=x, method=method, ...)
	}
)


# flag if one is executing gedProportions
.in_gedProportion <- sVariable(FALSE)

#' Estimating Cell Proportions from Known Signatures
#' 
#' \code{gedProportions} implements a pre-processing pipeline for applying deconvolution 
#' methodologies that use a known set of cell type-specific signatures in order to estimate 
#' cell proportions in heterogeneous samples (e.g., \cite{Abbas2009} or \cite{Gong2011}).
#' 
#' The actual estimation is performed via the \code{\link{ged}} interface, using a suitable 
#' deconvolution method. 
#' 
#' Before calling \code{\link{ged}}, the following pre-processing pipeline is applied to 
#' the data and/or the signature matrix:
#' \describe{
#' \item{map}{ the gene identifiers of the signature matrix into identifiers 
#' in the target global expression matrix;}
#' \item{subset}{ signatures and data matrices to a common set of features;}
#' \item{transform}{ signatures and data matrices to a common scale: linear or log;
#' Log-scale is automatically detected using the same heuristic as \code{GEO2R}.}
#' \item{normalise}{ jointly the signatures and data matrices using quantile normalisation.}  
#' }
#' 
#' All steps are optional and can be disabled if needed (see argument details).  
#' 
#' @param object target data, specified in any format supported by \code{\link{ged}}.
#' @param x basis signature data. 
#' It can be an \code{\linkS4class{ExpressionSet}} object or a \code{matrix}, whose columns 
#' contains cell-specific expression for each feature in the target data.
#'  
#' If the gene identifier type from the basis matrix do not match the one from the target matrix, 
#' these are converted using \code{\link{convertIDs}}.
#' If needed, this automatic conversion can be disabled using \code{map.method=NA}, as it is
#' by default when \code{x} is a \code{matrix}, whose rows are assumed to match the rows in the 
#' target matrix.
#' @param method method to use to deconvolve the target data and estimate cell proportions.
#' The method must be a deconvolution algorithm that is able to run using signatures as only auxiliary 
#' input.
#' 
#' The default method is \sQuote{lsfit}, which implements the algorithm proposed by \cite{Abbas2009}
#' that is based on standard regression.
#' An alternative method is the quadratic programming approach from \cite{Gong2011}, 
#' which solves a nonnegative least-square problem with sum-up-to one constraints on the 
#' proportions.
#' @param subset optional subset of features to use in the estimation. 
#' @param map.method method used to convert the basis signature's identifiers to match 
#' the target data's own type of identifiers. See \code{\link{mapIDs}}.
#' Identifier conversion can be disabled using \code{map.method=NA}. 
#' @param ... extra arguments passed to \code{\link{ged}}
#' @param log logical that indicates if the computation should take place in log or 
#' linear scale.
#' If \code{TRUE}, all non-log-scaled data (signatures and/or target) are log-transformed using
#' with base \code{lbase}.
#' If \code{FALSE}, all log-scaled data (signatures and/or target) are exp-transformed using
#' with base \code{lbase}.
#' If a number, then the function acts as if \code{log=TRUE} using the value of \code{log} 
#' as \code{lbase}. 
#' If \code{NULL}, then log-transform is applied only if either the signatures or the target 
#' data is in log scale, otherwise non-log-scaled data is exp-transformed into linear values, 
#' via \code{\link{expb}(A, lbase)}.
#' If \code{log=NA} no transformation is performed at all.
#' @param lbase numeric base use for the logarithmic/exponential transformations that 
#' are applied to the signature or data matrix.
#' @param normalize character string that specifies the normalisation method to apply 
#' jointly to the combined data (signatures,data).
#' The normalisation is performed after transforming the data and/or signatures if necessary.
#' @inheritParams ged
#' 
#' @seealso \code{\link{ged}}, \code{\link{gedBlood}} 
#' 
#' @importFrom preprocessCore normalize.quantiles
#' @export
gedProportions <- function(object, x, method='lsfit', subset=NULL, map.method=if( is.matrix(x) ) 'names' else 'auto'
						, ...
						, log=NULL, lbase=2, normalize=c('none', 'quantiles'), verbose=FALSE){
	
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	verbose <- lverbose()
	
	.in_gedProportion(TRUE)
	on.exit(.in_gedProportion(FALSE), add=TRUE)
	
	# match arguments
	if( isFALSE(normalize) || is_NA(normalize) ) normalize <- 'none'
	else if( isTRUE(normalize) ) normalize <- 'quantile'
	normalize <- match.arg(normalize)
	
	## WRAP INTO EXPRESSION SET(s)
	# signatures
	if( !is.matrix(x) && !is(x, 'ExpressionSet') )
		stop('Invalid basis signature object `x`: should be a matrix or an ExpressionSet object [', class(x), ']')
	doConvert <- TRUE
	if( is.matrix(x) ){
		idt <- asGeneIdentifierType(x)
		doConvert <- doConvert && !is.null(rownames(x)) && annotation(idt) != ''
#		vmessage("Wrapping signatures in an ExpressionSet object ... ", appendLF=FALSE)
#		x <- ExpressionSet(x, annotation=getAnnotation(x, null=FALSE))
#		vmessage('OK')
	}
	# data
	if( !is.matrix(object) && !is(object, 'ExpressionSet') )
		stop('Invalid target `object`: should be a matrix or an ExpressionSet object [', class(object), ']')
	if( is.matrix(object) ){
		idt <- asGeneIdentifierType(object)
		doConvert <- doConvert && !is.null(rownames(object)) && annotation(idt) != ''
#		vmessage("Wrapping target data in an ExpressionSet object ... ", appendLF=FALSE)
#		object <- ExpressionSet(object, annotation=getAnnotation(object, null=FALSE))
#		vmessage('OK')
	}
	##
	
	## SHOW initial features
	if( verbose > 1 ){
		vmessage("Original features:")
		vmessage(" * Data: ", str_out(featureNames(object)))
		vmessage(" * Signatures: ", str_out(featureNames(x)))
	}
	##
	sig <- x
	target <- object
	
	## ID CONVERSION
	vmessage("Mapping signature ids onto target ids ", appendLF=FALSE)
	if( doConvert && !is_NA(map.method) ){
		smap <- idFilter(map.method, .wrap=TRUE)
		vmessage("(method: ", smap$name, ") ... ", appendLF=FALSE)
		sig <- convertIDs(sig, target, method=map.method, verbose=verbose-1, rm.duplicates=TRUE)
		vmessage('OK [', nrow(sig), ' features x ', ncol(sig), ' cell types]')		
	}else{
		vmessage("(method: ", map.method, ") ... ", appendLF=FALSE)
		vmessage('... SKIP')
	}
	##
	
	## SUBSET: both the target and signature matrices
	vmessage("Limit/reorder to common set of features ... ", appendLF=FALSE)
	pids <- featureNames(sig)
	if( !is.null(pids) && !is.null(featureNames(target)) ){
		pids <- pids[pids %in% featureNames(target)]
		sig <- sig[pids, , drop=FALSE]
		target <- target[pids, , drop=FALSE]
		vmessage('OK [', nrow(sig), ' features x ', ncol(sig), ' cell types]')
	}else{
		vmessage('SKIP [no feature names]')
	}
	##
	
	## CHECK DIMS
	vmessage("Checking data dimension compatibility ... ", appendLF = FALSE)
	if( nrow(sig) != nrow(target) ){
		vmessage('ERROR')
		stop("Invalid basis signature matrix: number of rows [", nrow(sig), "] must be the same as in the data [", nrow(target), "]")
	}
	vmessage('OK [', nrow(sig), ' features x ', ncol(sig), ' cell types]')
	##
	
	## SUBSET-RESTRICT
	if( !is.atomic(subset) ){
    clsubset <- class(subset)
	  fsub <- featureNames(subset)
    if( hasAnnotation(subset) ){
      vmessage("Mapping subset ids onto target ids ", appendLF=FALSE)
      if( doConvert && !is_NA(map.method) ){
        smap <- idFilter(map.method, .wrap=TRUE)
        vmessage("(method: ", smap$name, ") ... ", appendLF=FALSE)
        subset <- convertIDs(subset, target, method=map.method, verbose=verbose-1)
        subset <- featureNames(subset)
        vmessage('OK [', length(fsub), ' -> ', length(subset), ']')
      }else{
        vmessage("(method: ", map.method, ") ... ", appendLF=FALSE)
        vmessage('... SKIP')
        subset <- featureNames(subset)
      }
    }
		if( is.null(subset) )
			warning("Not using subset specification: extracted NULL feature names from object of class '", clsubset, "'")
	}
	
	if( !is.null(subset) ){
		vmessage("Limiting features to subset (type: ", class(subset), " [", str_out(subset), "]) ... ", appendLF=FALSE)
		# do subset
		if( is.character(subset) ){
			target <- target[featureNames(target) %in% subset, , drop=FALSE]
			sig <- sig[featureNames(sig) %in% subset, , drop=FALSE]
		}else{
			target <- target[subset, , drop=FALSE]
			sig <- sig[subset, , drop=FALSE]
		}
		vmessage('OK [', nrow(sig), ' features x ', ncol(sig), ' cell types]')
		
		## CHECK DIMS
		vmessage("Re-checking dimension compatibility ... ", appendLF = FALSE)
		if( nrow(sig) == 0 || nrow(target) == 0 ){
		  vmessage('ERROR')
		  stop("Empty data after subsetting features: signatures [", nrow(sig), "] - target data [", nrow(target), "]")
		}
    
		if( nrow(sig) != nrow(target) ){
			vmessage('ERROR')
			stop("Invalid basis signature matrix after subsetting: number of rows [", nrow(sig), "] must be the same as in the data [", nrow(target), "]")
		}
		vmessage('OK')
		##
	}
	##
	
	## SHOW feature names
	if( verbose > 1 ){
		vmessage("Using features:")
		vmessage(" * Data: ", str_out(featureNames(target)))
		vmessage(" * Signatures: ", str_out(featureNames(sig)))
	}
	##
	
	## TRANSFORM: log-transform the signatures if necessary
	# data
	vmessage("Checking log-scale ... ", appendLF = FALSE)
	log_data <- is_logscale(target)
	vmessage("data:", if( log_data ) "YES" else "NO", appendLF=FALSE)
	# signatures
	log_sig <- is_logscale(sig)
	vmessage(" - signatures:", if( log_sig ) "YES" else "NO")
	
	if( is.null(log) ){ # log transform if either the signatures or the data are in log scale
		log <- log_sig || log_data
	}else if( isNumber(log) ){
		lbase <- log
		log <- TRUE
	}
	
	if( isTRUE(log) ){ # LOG
		if( !log_sig ){
			vmessage("Applying log-transform to signatures (base ", lbase, ") ... ", appendLF = FALSE)
			sig <- log_transform(sig, lbase)
			vmessage("OK")
		}	
		if( !log_data ){
			vmessage("Applying log-transform to data (base ", lbase, ") ... ", appendLF = FALSE)
			target <- log_transform(target, lbase)
			vmessage("OK")
		}
	} else if( !is_NA(log) ){ # LINEAR
		if( log_sig ){
			vmessage("Reverting log-transform on signatures (base ", lbase, ") ... ", appendLF = FALSE)
			sig <- expb(sig, lbase)
			vmessage("OK")
		}	
		if( log_data ){
			vmessage("Reverting log-transform on data (base ", lbase, ") ... ", appendLF = FALSE)
			target <- expb(target, lbase)
			vmessage("OK")
		}
	}
	##
	
	# NORMALISE: normalise both data together if necessary
	vmessage("Normalizing signatures and target together (method: ", normalize, ") ... ", appendLF = FALSE)
	if( normalize != 'none' ){
		mSig <- exprs(sig)
		## combine data
		xAll <- cbind(mSig, exprs(target))
		allNorm <- 
		switch(normalize,
			quantiles = preprocessCore::normalize.quantiles(xAll)
			, stop("Unexpected error: unsupported normalization method '", normalize, "'")
		)
		vmessage("OK")
		
		## split normalized data (restoring dimnames)
		# signatures
		sigN <- allNorm[, 1:ncol(mSig), drop=FALSE]
		dimnames(sigN) <- dimnames(mSig)
		if( is(sig, 'ExpressionSet') ){
			exprs(sig) <- sigN
			validObject(sig)
		}else sig <- sigN
		# target
		xN <- allNorm[, (ncol(sig)+1):ncol(allNorm), drop=FALSE]
		if( is(target, 'ExpressionSet') ){
			dimnames(xN) <- dimnames(exprs(target))
			exprs(target) <- xN
		}else{
			dimnames(xN) <- dimnames(target)
			target <- xN
		}
		validObject(target)
	}else{
		vmessage("SKIP")
	}
	
	# GED: apply ged method
	res <- ged(target, sig, method=method, ..., verbose = verbose)
	
	# return res
	res
	
}

#' Checking a Deconvolution Method
#' 
#' This function runs a given algorithm on a small random toy dataset, 
#' to check that the algorithm works correctly
#' 
#' @param method algorithm specification
#' @param ... other arguments passed to \code{\link{ged}}.
#' @param maxIter maximum number of iterations.
#' Note that it is fixed to a small default value for the purpose of the check, 
#' which is likely to make no sense for real application.
#' @param verbose defines verbosity level (\code{FALSE/TRUE} or a number).
#' @return returns invisibly the result of running the algorithm, which 
#' usually is an \code{NMFfit} object.
#' 
#' @export
#' @examples
#' gedCheck('deconf')
#' gedCheck('DSA', log=FALSE)
#' 
gedCheck <- function(method, ..., maxIter=5L, verbose=3L){
	
	# random data
	x <- rmix(3, 100, 20)
	model <- mixData(x)
	marks <- getMarkers(x)
	
	# load algorithm
	meth <- gedAlgorithm(method)
	
	# apply algorithm on relevant input data
	res <- if( onlyRequired('Basis', meth) ) ged(x, basis(model), meth, maxIter=maxIter, ..., verbose=verbose)
			else if( onlyRequired('Coef', meth) ) ged(x, coef(model), meth, maxIter=maxIter, ..., verbose=verbose)
			else if( onlyRequired('Marker', meth) ) ged(x, marks, meth, maxIter=maxIter, ..., verbose=verbose)
			else if( !anyRequired(meth) ) ged(x, nbasis(model), meth, maxIter=maxIter, ..., verbose=verbose)
	
	if( !isNMFfit(res) ) message('NOTE: result of ged was not an NMFfit object [', class(res), ']') 
	invisible(res)
}
