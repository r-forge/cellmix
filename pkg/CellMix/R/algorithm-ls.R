# Deconvolution algorithms based on standard least-squares:
#
# Simple strategy from Abbas et al. (2009)
# 
# Author: Renaud Gaujoux
# Creation: 18 Jan 2012
###############################################################################

#' @include registry-algorithms.R
NULL

# lsfit + iterations to enforce nonnegative values
.nn_lsfit <- function(x, y, nneg=c('iterate', 'pmax', 'none'), ...){
	
	if( !is_NA(nneg) )
		nneg <- match.arg(nneg)
	
	# fit with lsfit
	fit <- lsfit(x=x, y=y, ...)
	# store returned values from lsfit
	res <- list()
	res$lsfit <- fit
	
	# enforce nonnegative values if necessary
	res$coef <- fit$coef
	if( !is_NA(nneg) ){
		if( nneg == 'pmax' ) res$coef <- pmax(fit$coef) 
		else if( nneg == 'iterate' ){
			r <- ncol(x)
			# iterate as described in Abbas et al. (2009)
			res$coef <- 
			sapply(seq(ncol(fit$coef)), function(k){
				a <- fit$coef[,k]
				i <- integer()
				while( any(a<0) && length(i) < r ){
					i <- c(i, which.min(a))
					a[i] <- 0
					fit <- lsfit(x=x[, -i, drop=FALSE], y=y[,k], ...)
					stopifnot( length(a[-i]) == length(fit$coef) )
					a[-i] <- fit$coef
				}
				a
			})
		}
	}
	
	# re-apply dimnames
	rownames(res$coef) <- colnames(x)
	colnames(res$coef) <- colnames(y)
	
	res
}

#' Partial Gene Expression Deconvolution by Least-Square 
#' 
#' @inheritParams .qprog
#' @param nneg specification of the method used to enforce the nonnegativity 
#' of the estimated proportions (used only in method \dQuote{lsfit}).
#' @param rescale logical used when estimating proportions from signatures, 
#' that indicates if the esti,ated coefficients should be scaled to sum up to 
#' one (\code{TRUE}) or left as estimated by the linear regression (\code{FALSE}).
#' This scaling is performed after the coefficients have been forced to be 
#' nonnegative.
#' @param ... extra arguments passed to \code{\link{lsfit}}.   
#' 
#' @return an \code{\link{NMF}} object.
#' 
#' @keywords internal
.gedLSfit <- function(X, seed, nneg=c('iterate', 'pmax', 'none'), rescale=TRUE, ...){
	
	if( hasBasis(seed) && !hasCoef(seed) ){
		
		nneg <- match.arg(nneg)
		vmessage("Estimating cell proportions from cell-specific signatures [lsfit]")
		# fit with lsfit
		res <- .nn_lsfit(x=basis(seed), y=X, intercept=FALSE, ..., nneg=nneg)
		# store returned values from lsfit
		seed$lsfit <- res$lsfit		
		# update coef
		coef(seed) <- res$coef
		# rescale coef
		if( rescale ) coef(seed) <- scoef(seed)
		
	}else if( !hasBasis(seed) && hasCoef(seed) ){
		
		vmessage("Estimating cell-specific signatures from cell proportions [cs-lsfit]")
		
		# by default enforce nonnegativity only if the input data contains negative values 
		if( missing(nneg) && min(X) < 0 ) nneg='none'

		# fit the transposed problem with lsfit
		res <- .nn_lsfit(x=t(coef(seed)), y=t(X), intercept=FALSE, ..., nneg=nneg)
		# store returned values from lsfit
		seed$lsfit <- res$lsfit
		
		# update basis 
		basis(seed) <- t(res$coef)
		
	}else{
		stop("lsfit - Only partial fitting is currently supported (i.e. either the signatures or the proportions must be provided).")
	}
	
	# return updated seed
	seed
}

# Registration of NMF method 'lsfit'
nmfAlgorithm.lsfit <- setNMFMethod('lsfit', .gedLSfit,
	objective='euclidean', 
	mixed=TRUE,
	overwrite=TRUE
)

#' Partial Gene Expression Deconvolution by Standard Least-Squares
#' 
#' Estimates cell/tissue proportions given a known set of cell/tissue-specific 
#' expression signatures, using standard least-squares as proposed by 
#' \cite{Abbas2009}.
#'
#' The algorithm uses a heuristic to enforce the nonnegativity of the estimated 
#' proportions, that consists in fitting successive regressions, each time excluding 
#' the most negative coefficient from the model, until all coefficients are nonnegative.
#'
#' All regressions are fitted using the function \code{\link{lsfit}}.
#'  
#' @inheritParams .gedLSfit
#' @aliases lsfit-ged
#' @cite Abbas2009
#' @examples 
#' 
#' # random target matrix
#' x <- rmatrix(100, 20)
#' # random cell signatures
#' s <- rmatrix(100, 3)
#' 
#' # deconvolve using standard least-squares
#' res <- ged(x, s, 'lsfit')
#' coef(res)
#' # signatures are not updated
#' identical(basis(res), s)
#' \dontshow{ 
#' 	stopifnot(identical(basis(res), s))
#'	stopifnot( nmf.equal(res, ged(x, s, 'lsfit')) ) 
#' }
#' 
#'
gedAlgorithm.lsfit <- setGEDMethod(key='lsfit'
	, description = "Partial deconvolution of proportions using standard least-squares"
	, algorithm = 'lsfit' 
	, reqBasis = TRUE
	, reqCoef= FALSE
	, reqMarker = FALSE
	, maxIter = 1L
	, cite = "Abbas2009"
)

#' Cell-Specific Expression by Standard Least-Squares
#' 
#' Estimates cell-specific proportions given known proportions 
#' expression signatures, using standard least-squares.
#'
#' The algorithm uses the same heuristic as the ged algorithm \code{\link[=lsfit-ged]{lsfit}} 
#' to enforce non-negative expression values.
#' It is included in the \pkg{CellMix} package for test/experimental purposes. 
#'  
#' @inheritParams gedAlgorithm.lsfit
#' @aliases cs-lsfit-ged
#' @examples 
#' 
#' # random target matrix
#' x <- rmatrix(100, 20)
#' # random cell proprtions
#' p <- rmatrix(3, 20)
#' 
#' # deconvolve using standard least-squares
#' res <- ged(x, p, 'cs-lsfit')
#' head(basis(res))
#' # proportions are not updated
#' identical(coef(res), p)
#' \dontshow{ 
#' 	stopifnot(identical(coef(res), p))
#'	stopifnot( nmf.equal(res, ged(x, p, 'cs-lsfit')) ) 
#' }
#' 
gedAlgorithm.cs_lsfit <- setGEDMethod(key='cs-lsfit'
		, description = "Partial deconvolution of cell signatures using standard least-squares"
		, algorithm = 'lsfit'
		, reqBasis = FALSE
		, reqCoef= TRUE
		, reqMarker = FALSE
		, maxIter = 1L
)

#' Partial Gene Expression Deconvolution by Nonnegative Least-Square 
#' 
#' @inheritParams .gedLSfit
#' @param ... extra arguments passed to \code{\link[NMF]{.fcnnls}}.
#' @return an \code{\link{NMF}} object.
#' 
#' @keywords internal
.gedNNLS <- function(X, seed, rescale=TRUE, ...){
	
	if( hasBasis(seed) && !hasCoef(seed) ){
		
		vmessage("Estimating cell proportions from cell-specific signatures [fcnnls]")
		# fit with lsfit
		res <- .fcnnls(basis(seed), X, ...)
		# update coef
		coef(seed) <- res$coef
		# rescale coef
		if( rescale ) coef(seed) <- scoef(seed)
		
	}else if( !hasBasis(seed) && hasCoef(seed) ){
		
		vmessage("Estimating cell-specific signatures from cell proportions [fcnnls]")
		# fit the transposed problem with lsfit
		res <- .fcnnls(t(coef(seed)), t(X), ...)
		# update basis
		basis(seed) <- t(res$coef)
		
	}else{
		stop("fcnnls - Only partial fitting is currently supported (i.e. either the signatures or the proportions must be provided).")
	}
	
	# return updated seed
	seed
}

# Registration of NMF method 'lsfit'
setNMFMethod('nnls', .gedNNLS,
		objective='euclidean',
		overwrite=TRUE
)


#' Partial Gene Expression Deconvolution by Nonnegative Least-Squares
#' 
#' Estimates cell/tissue proportions given a known set of cell/tissue-specific 
#' expression signatures, using the fast combinatorial nonnegative least-square 
#' method, which was proposed by \cite{KimH2007} to perform nonnegative matrix 
#' factorization of gene expression -- not originally for deconvolution.
#'
#' All regressions are fitted using the function \code{\link[NMF]{.fcnnls}}.
#'  
#' @inheritParams .gedNNLS
#' @aliases nnls-ged
#' @examples 
#' 
#' # random target matrix
#' x <- rmatrix(100, 20)
#' # random cell signatures
#' s <- rmatrix(100, 3)
#' 
#' # deconvolve using nonnegative least-squares
#' res <- ged(x, s, 'nnls')
#' coef(res)
#' # signatures are not updated
#' identical(basis(res), s)
#' \dontshow{ 
#' 	stopifnot(identical(basis(res), s))
#'	stopifnot( nmf.equal(res, ged(x, s, 'nnls')) ) 
#' }
#' 
#'
gedAlgorithm.nnls <- setGEDMethod(key='nnls'
		, description = "Partial deconvolution of proportions or signatures using fast combinatorial nonnegative least-squares"
		, algorithm = 'nnls' 
		, reqBasis = TRUE
		, reqCoef= FALSE
		, reqMarker = FALSE
		, maxIter = 1L
)

#' Partial Cell-Specific Expression by Nonnegative Least-Squares
#' 
#' Estimates cell-specific signatures given known cell proportions, 
#' using the fast combinatorial nonnegative least-square as the ged algorithm 
#' \code{\link[=nnls-ged]{nnls}}.
#' It is included in the \pkg{CellMix} package for test/experimental purposes.
#'
#' All regressions are fitted using the function \code{\link[NMF]{.fcnnls}}.
#'  
#' @inheritParams .gedNNLS
#' @aliases cs-nnls-ged
#' @examples 
#' 
#' # random target matrix
#' x <- rmatrix(100, 20)
#' # random cell proprtions
#' p <- rmatrix(3, 20)
#' 
#' # deconvolve using nonnegative least-squares
#' res <- ged(x, p, 'cs-nnls')
#' head(basis(res))
#' # proportions are not updated
#' identical(coef(res), p)
#' \dontshow{ 
#' 	stopifnot(identical(coef(res), p))
#'	stopifnot( nmf.equal(res, ged(x, p, 'cs-nnls')) ) 
#' }
#' 
gedAlgorithm.cs_nnls <- setGEDMethod(key='cs-nnls'
		, description = "Partial deconvolution of signatures using fast combinatorial nonnegative least-squares"
		, algorithm = 'nnls'
		, reqBasis = FALSE
		, reqCoef= TRUE
		, reqMarker = FALSE
		, maxIter = 1L
)
