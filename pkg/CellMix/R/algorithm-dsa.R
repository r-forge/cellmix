# Digital Sorting Algorithm
# 
# Author: Renaud Gaujoux
# Created: Mar 26, 2013
###############################################################################

#' @include algorithm-qprog.R
NULL

.data2exp <- function(X, log = NULL){
	vmessage('Checking data scale ... ', appendLF=FALSE)
	
	if( !is_logscale(X) ){
		vmessage('OK [linear]')
		return(X)
	}else{
		vmessage('NOTE [log]')
		vmessage('Converting data to linear scale ... ', appendLF=FALSE)
		
		# skip if requested
		if( isFALSE(log) ){
			vmessage('SKIP')
			return(X)
		}
		
		lbase <- 2
		if( isNumber(log) ){
			lbase <- log
		}
		
		# convertion
		X <- expb(X, lbase)
		vmessage('OK [base: ', lbase,']')
	}
	# converted data
	X
}

#' Digital Sorting Algorithm: Proportion Estimation Method   
#' 
#' \code{DSAproportions} implements the method proposed by 
#' \cite{Zhong2013} as part of the \code{\link[=DSA-ged]{DSA}} 
#' algorithm.
#' This method estimates cell proportions from mixed sample expression data, 
#' given a set of markers, i.e. features that are known to be 
#' exclusively expressed by a single cell type in the mixture.  
#' 
#' @param x target mixed expression matrix
#' @param markers marker specification
#' @param log indicates if the data are in log-scale or should be converted to linear-scale.
#' This is relevant because the DSA algorithm assumes that the input mixed data are in linear 
#' scale (i.e. not log-trasnformed).
#' 
#' If \code{NULL}, then data's scale is detected by \code{link{is_logscale}} 
#' and conversion to linear-scale is performed if necessary.
#' If \code{TRUE} the data is exponentialised (using log base 2).
#' If \code{FALSE} the data is left unchanged (the detected log scale is
#' displayed in verbose mode).
#' If a number, then it is used as the log base to exponentialise the data.
#' @param match.names logical that indicates if marker names should be first matched
#' against the data, to convert them into integer indexes.
#' If code{FALSE}, then \code{markers} should be a list of integer vectors that 
#' corresponds to the indexes of markers for each cell types.
#' @param verbose verbose level
#' 
#' @return a matrix of proportions with the same number of 
#' columns as \code{x} and as many rows as elements  in 
#' \code{markers} (i.e. cell types).
#' 
#' @export
#' @examples
#' 
#' # random data
#' x <- rmix(3, 100, 20)
#' m <- getMarkers(x)
#' 
#' # estimate proportions
#' p <- DSAproportions(x, m)
#' # plot against true proportions
#' profplot(x, p)
#' 
DSAproportions <- function(x, markers, log = NULL, match.names = TRUE, verbose=FALSE){
	
	# set verbosity
	ov <- lverbose(verbose)
	on.exit( lverbose(ov) )
	
	if( !is.list(markers) ){
		stop("Invalid argument 'markers': must be a list of markers.")
	}
	idx <- if( match.names ) matchIndex(markers, x) else markers
	
	# IMPORTANT: convert to linear scale
	x <- .data2exp(x, log = log)
	
	# compute proportions
	vmessage("Computing proportions using DSA method ... ", appendLF=FALSE)
	p <- .DSAproportions(exprs(x), idx)
	vmessage("OK")
	
	p
}
#' \code{.DSAproportions} is the internal -- non-exported -- 
#' function that actually performs the estimation.
#' It expects a list of valid indexes for each cell type
#' whose proportions are to be estimated. 
#'  
#' @rdname DSAproportions 
.DSAproportions <- function(x, markers){
		
	# compute average expression
	O <- colMeansBy(x, markers)
	# solve linear equations
#	f <- lsfit(t(O), rep(1, ncol(O)), intercept=FALSE)
#	g_ <- coef(f)
	f <- fcnnls(t(O), matrix(1, ncol(O), 1))
	g_ <- as.numeric(f$x)
	p <- diag(g_) %*% O
	# add rownames if necessary
	if( is.null(rownames(p)) ){
		rownames(p) <- rownames(O)
	}
	
	# return proportion estimate
	p
}

#' Complete Deconvolution using Digital Sorting Algorithm (DSA) 
#' 
#' The method \dQuote{DSA} implements the Digital Sorting Algorithm (DSA) 
#' proposed by \cite{Zhong2013}, which performs complete gene 
#' expression deconvolution using a set of marker genes only.
#'  
#' @inheritParams .qprog
#' @aliases DSA-ged
#' @examples 
#' 
#' ## Example on dummy/random data
#' x <- rmix(3, 100, 20)
#' # markers
#' ml <- getMarkers(x)
#' 
#' # deconvolve using DSA (quadratic programming)
#' res <- ged(x, ml, 'DSA', verbose=TRUE, log=FALSE)
#' # plot against true data
#' profplot(x, res)
#' profplot(t(basis(x)), t(basis(res)), legend=FALSE)
#' 
#' \dontshow{ 
#' 	stopifnot( nmf.equal(res, ged(x, ml, 'DSA', log=FALSE), identical=FALSE) ) 
#' }
#' 
#' @demo DSA Algorithm: Digital Sorting Algorithm
#' ####################################################
#' # Sample deconvolution analysis with DSA in CellMix  
#' ####################################################
#' 
#' # load benchmark data
#' x <- ExpressionMix('GSE19830', verbose=TRUE)
#' dim(x)
#' annotation(x)
#' # extract mixed samples
#' mix <- mixedSamples(x)
#' 
#' # load TIGER marker list
#' ml <- MarkerList('TIGER')
#' ml
#' names(ml)
#' 
#' # select markers for the tissues present in the mixture
#' basisnames(x)
#' ml <- ml[c('brain', 'liver', 'lung')]
#' summary(ml)
#' 
#' # convert to match annotations
#' mlx <- convertIDs(ml, mix, verbose=TRUE)
#' summary(mlx)
#' 
#' # QC on markers from their expression patterns in mixed samples
#' profplot(mlx[,1:10], mix)
#' # filter out poor markers using SCOREM (based on linear-scale expression)
#' mlsc <- extractMarkers(mlx, expb(mix, 2), method='SCOREM', alpha=10^-12)
#' summary(mlsc)
#' # expresison patterns are more correlated
#' profplot(mlsc[,1:10], mix)
#' 
#' # apply DSA using all markers
#' res <- ged(mix[mlsc,], mlsc, 'DSA', verbose=TRUE)
#' 
#' # plot against true proportions
#' profplot(mix, res)
#' 
gedAlgorithm.DSA <- setGEDMethod(key='DSA'
		, description = "Complete deconvolution using Digital Sorting Algorithm"
		, algorithm = 'qprog'
		, reqBasis = FALSE
		, reqCoef= FALSE
		, reqMarker = TRUE
		, maxIter = 1L
		, cite = 'Zhong2013'
		, defaults = list(seed='none')
)
