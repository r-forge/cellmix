# Digital Sorting Algorithm
# 
# Author: Renaud Gaujoux
# Created: Mar 26, 2013
###############################################################################

#' @include algorithm-qprog.R
NULL


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
DSAproportions <- function(x, markers, verbose=FALSE){
	
	if( verbose ) message("Computing proportions using DSA method ... ", appendLF=FALSE)
	
	if( !is.list(markers) ){
		stop("Invalid argument 'markers': must be a list of markers.")
	}
	idx <- matchIndex(markers, x)
	p <- .DSAproportions(exprs(x), idx)
	
	if( verbose ) message("OK")
	
	p
}
#' \code{.DSAproportions} is the internal -- non-exported -- 
#' function that actually performs the estimation.
#' It expects a list of valid indexes for each cell type
#' whose proportions are to be estimated. 
#'  
#' @rdname DSAproportions 
.DSAproportions <- function(x, markers, verbose=FALSE){
	
	if( verbose ) message("Computing proportions using DSA method ... ", appendLF=FALSE)
	
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
	
	if( verbose ) message("OK")
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
