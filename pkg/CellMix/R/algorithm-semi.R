# Semi-supervised algorithms for gene expression deconvolution
# 
# Author: Renaud Gaujoux
# Creation: 22 Jan 2012
###############################################################################

#' @include algorithm-common.R
#' @include simulations.R
NULL


# Expose non exported functions and methods	
name <- NMF:::name


###############################################################################
# Define new/modified NMF algorithms that use the marker list to ensure that
# the components are identifiable.
###############################################################################

#' Initialisation Method for semi-NMF Algorithms
#' 
#' Pre-process and attach markers to an NMF model so that they can be enforced at
#' each iteration.
#' It adds a static variable 'markerIdx' to the workspace.
#' 
#' @param method method name, used in case of error only
#' @param v target matrix
#' @param x initial NMF model 
#' @param data marker list
#' @param ratio expression ratio of markers between its cell type and 
#' other cell types.
#' @param eps lower limit for cell type signature expression levels or proportions.
#' This is used to avoid null values that might be problematic for some deconvolution algorithm, 
#' e.g., for the method \sQuote{ssKL} that is based on the KL-divergence, which computes logs.
#' @inheritParams .gedlab
#' 
#' @return the modified NMF model. 
#' 
#' @keywords internal
markerInit <- function(method, v, x, data=NULL
						, markers=c('prior+semi', 'semi')
						, ratio=NULL, eps=NULL, ...){

	nm <- if( !is.null(method) ) name(method)
	verbose <- nmf.getOption('verbose')
	
	if( verbose ) message("Checking for markers in call ... ", appendLF=FALSE)
	if( is.null(data) ){
		staticVar('markerIdx', NULL)
		if( verbose ) message("NO [unsupervised]")
		return(x)
	}
	markers <- match.arg(markers)
	if( verbose ) message("OK [", markers, ']')
	
	# cs-fold
	csfold <- if( is.null(ratio) ) Inf else ratio
	if( verbose ) message("Using cell-specific marker expression fold: ", csfold)
	
	# get the marker indexes
	if( verbose ) message("Matching markers against data ... ", appendLF=FALSE)
	M <- matchIndex(data, v)
	if( !all(sapply(M, is.integer)) ){
		if( verbose ) message("ERROR")
		stop('Unexpected error - invalid marker list: matched marker IDs must be integer indexes.')
	}
	if( verbose ) message("OK [", nmark(M), "/", nmark(data), ']')
	
	# define how to use markers 
	doSeed <- grepl('prior', markers)
	
	# verbose message
	nmN <- nmark(M, each=TRUE)
	nml <- sum(nmN)
	if( !nml ){
		stop("cannot apply ged method '", nm, "': could not match any markers [", str_out(marknames(data)), '] in target data [', str_out(featureNames(v)), ']')
	}
	if( any(w <- nmN == 0) ){
		stop("cannot apply ged method '", nm, "': could not match any markers for types ", str_out(names(data)[w]), ' in target data')
	}
	
	# seed proportions using method from DSA [Zhong et al. (2013)]
	if( doSeed || !hasCoef(x) ){
		p <- .DSAproportions(v, M, verbose=verbose)
		if( !is.null(eps) ){ # shunt small values
			p <- pmax(p, eps)
		}
		coef(x) <- p
	}
	if( doSeed || !hasBasis(x) ){
		csfold1 <- if( is.null(ratio) ) 1.5 else ratio
		basis(x) <- csfit_markers(v, coef(x), M, ratio=csfold1, verbose=verbose)		
	}
	# enforce an initial block structure on the model
	if( verbose ) message("Enforcing marker initial expresison pattern (cs-fold: ", csfold, ") ... ", appendLF=FALSE)
	x <- enforceMarkers(x, M, ratio=ratio, value=max(basis(x)), eps=eps, attach=TRUE)
	if( verbose ) message("OK")
	
	# store the markers in a static variable for possible enforcing at each iteration
	staticVar('markerIdx', M)
	
	# return modified model
	x
}

###################################################
# Define modified version of Brunet's algorithm
###################################################
nmfAlgorithm.ssKL <- setNMFMethod('ssKL', 'KL'
	, onInit = markerInit
	, onReturn = function(x){ 
		coef(x) <- scoef(x)
		x
	}
	, defaults = list(seed='rprop', eps=.Machine$double.eps)
	, Update=function(i, v, x, copy=FALSE, data=NULL, sscale=FALSE, ratio=NULL, alpha=0, ...)
	{	
		# setup marker index vector
		if( i == 1 ){			
			# enforce the copy for the first iteration
			#copy <- TRUE
			#x <- markerInit('Gbrunet', v, x, data, ratio)
			
			# add soft constraint on sum to one
			if( alpha > 0 )
				staticVar('v', rbind(v, rep(alpha,ncol(v))))
		}
		
		# retrieve each factor
		w <- basis(x); h <- coef(x);
		if( alpha > 0 ){
			copy <- TRUE
			w <- rbind(w, rep(alpha, nbasis(x)))
			v <- staticVar('v')
		}
		
		# standard divergence-reducing NMF update for H	
		h <- nmf_update.KL.h(v, w, h, copy=copy)		
		
		# standard divergence-reducing NMF update for W
		w <- nmf_update.KL.w(v, w, h, copy=copy)
		
		#every 10 iterations: adjust small values to avoid underflow
		# NB: one adjusts in place even when copy=TRUE, as in this case 'h' and 'w' 
		# are mapped to locally allocated memory
		if( i %% 10 == 0 ){
			eps <- .Machine$double.eps
			# apply min for numerical stability
			h <- pmax.inplace(h, eps)
			w <- pmax.inplace(w, eps)
		}
		# apply in place inequality constrains from the markers
		if( !is.null(data) )
			w <- neq.constraints.inplace(w, staticVar('markerIdx'), ratio=ratio)
		
		# update object if the updates duplicated the data
		if( copy ){		
			#return the modified data	
			.basis(x) <- if( alpha > 0 ) w[-nrow(w),] else w; 
			.coef(x) <- h;
		}
		
		# scale basis components
		if( sscale ){
			a <- colMeans(.basis(x), na.rm=TRUE)
			.basis(x) <- sweep(.basis(x), 2L, a, "/", check.margin=FALSE)
			if( sscale == 1 )
				.coef(x) <- sweep(.coef(x), 1L, a, "*", check.margin=FALSE)
		}
		
		return(x)
	}
, overwrite=TRUE)


#' Complete Gene Expression Deconvolution by Semi-Supervised NMF
#' 
#' Algorithms \sQuote{ssKL} and \sQuote{ssFrobenius} are modified versions of the
#' original NMF algorithm from \cite{Brunet2004} and \cite{Lee2001}, 
#' that use a set of known marker genes for each cell type, to enforce the 
#' expected block expression pattern on the estimated signatures, as 
#' proposed in \cite{Gaujoux2011}.
#' 
#' These algorithms simultaneously estimates both the cell-specific signature and mixture 
#' proportion matrices, using block-descent method that alternately estimates each matrix.
#' Both re-scale the final proportion estimates so that they sum-up to one. 
#' 
#' The functions \code{gedAlgorithm.ssKL} and \code{gedAlgorithm.ssFrobenius} are wrapper
#' functions to the underlying NMF algorithms.
#' They are primiraly defined to enable correct listing their specific arguments on this 
#' help page. 
#' The recommend way of applying these algorithms is via \code{\link{ged}} interface 
#' (e.g., \code{ged(..., method='ssKL')}).
#' 
#' @inheritParams NMF::nmf.stop.stationary
#' @inheritParams NMF::nmfAlgorithm.KL
#' @inheritParams markerInit
#' @param seed default seeding method.
#' @param sscale specifies how signatures -- and proportions -- are re-scaled at the end 
#' of each iteration.
#' If \code{TRUE}, each signature is mean-centered separately.
#' If \code{2}, then each signature is mean-centered separately and the inverse linear transformation 
#' proportions is applied to the proportions (i.e. on the rows of the mixture coefficient matrix), so 
#' that the fitted matrix does not change.
#' If \code{FALSE}, no re-scaling is performed at all.
#' @param alpha numeric coefficient used to smoothly enforce a sum-up-to-one constraint on the 
#' proportions, by regularising the objective function.
#' If \code{NULL,} no constraint is applied.  
#' 
#' @aliases ssKL-ged
#' @cite Gaujoux2011
#' 
gedAlgorithm.ssKL <- setGEDMethod(key='ssKL'
	, description = "Semi-supervised NMF algorithm for KL divergence, using marker genes"
	, algorithm = 'ssKL'
	, reqBasis = FALSE
	, reqCoef= FALSE
	, reqMarker = TRUE
	, maxIter = 3000L
	, cite = "Gaujoux2011"
)

#######################################################
# Define modified version of Lee and Seung's algorithm
#######################################################
nmfAlgorithm.ssFrobenius <- setNMFMethod('ssFrobenius', 'Frobenius'
	, defaults = list(seed='rprop')
	, onInit = markerInit # seeding using markers 
	, onReturn = function(x){ 
		# scale to proportions
		coef(x) <- scoef(x)
		x
	}
	, Update = function(i, v, x, sscale=TRUE, copy=FALSE, data, ratio=NULL, alpha=0, ...)
	{
		# setup marker index vector
		if( i == 1 ){
			# add soft constraint on sum to one
			if( alpha > 0 )
				staticVar('v', rbind(v, rep(alpha, ncol(v))))
		}
		
		# retrieve each factor
		w <- basis(x); h <- coef(x);	
		if( alpha > 0 ){
			copy <- TRUE
			w <- rbind(w, rep(alpha, nbasis(x)))
			v <- staticVar('v')
		}
		
		#precision threshold for numerical stability
		eps <- 10^-9
		
		# euclidean-reducing NMF iterations	
		# H_au = H_au (W^T V)_au / (W^T W H)_au
		h <- nmf_update.euclidean.h(v, w, h, eps=eps, copy=copy)
		# update original object if not modified in place
		if( copy ) .coef(x) <- h
		
		# W_ia = W_ia (V H^T)_ia / (W H H^T)_ia and columns are rescaled after each iteration
		w <- nmf_update.euclidean.w(v, w, h, eps=eps, copy=copy)
		
		if( alpha > 0 ) w <- w[-nrow(w),]
		
		# apply in place inequality constrains from the markers
		if( !is.null(data) ){
			w <- neq.constraints.inplace(w, staticVar('markerIdx'), ratio=ratio)
		}
		
		#rescale columns TODO: effect of rescaling? the rescaling makes the update with offset fail
		if( sscale ){
			a <- colMeans(w, na.rm=TRUE)
			w <- sweep(w, 2L, a, "/", check.margin=FALSE)
			if( sscale == 2 )
				.coef(x) <- sweep(h, 1L, a, "*", check.margin=FALSE)
		}
		
		#return the modified data		
		.basis(x) <- w
		return(x)
	}
	, overwrite = TRUE
)

#' @inheritParams NMF::nmfAlgorithm.Frobenius
#' @rdname gedAlgorithm.ssKL
#' @aliases ssFrobenius-ged
gedAlgorithm.ssFrobenius <- setGEDMethod(key='ssFrobenius'
		, description = "Semi-supervised NMF algorithm for Euclidean distance, using marker genes"
		, algorithm = 'ssFrobenius'
		, reqBasis = FALSE
		, reqCoef= FALSE
		, reqMarker = TRUE
		, maxIter = 3000L
		, cite = "Gaujoux2011"
)
