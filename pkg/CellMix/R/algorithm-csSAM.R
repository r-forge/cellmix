# csSAM wrapper
# 
# Wraps call to the csSAM package from Shen-Orr et al. (2010)
# Author: Renaud Gaujoux
###############################################################################

#' @include registry-algorithms.R
#' @include plots.R
NULL

# simple warning note
wnote <- function(..., immediate. = TRUE){
	warning(..., immediate.=immediate., call.=FALSE)
}

#' Cell-specific Differential Expression with csSAM
#' 
#' This function is adapted from the function \code{\link[csSAM]{csSamWrapper}} 
#' in the \code{\link{csSAM}} package, to integrate the csSAM algorithm from 
#' \cite{Shen-Orr2010} into the \pkg{CellMix} framework of deconvolution algorithms.
#' 
#' @inheritParams csSAM::csSamWrapper
#' @param Y target global gene expression matrix (n x p), with samples in columns, ordered in the 
#' same order at the cell proportions data in \var{x}.
#' @param x known cell proportions as a matrix (k x p) or an \code{\linkS4class{NMF}} model 
#' containing the cell proportions in the coefficient matrix -- and a normally
#' empty basis matrix.
#' The proportions must be ordered in the same order as the samples in the target matrix.
#' 
#' For \code{csTopTable}, a csSAM fit as return by \code{\link{ged}}.
#' @param data specification of the sample groups.
#' If not missing, it must be a factor or coercible to a factor, with length the 
#' number of samples, i.e. columns, in the target matrix. 
#' @param nperms The number of permutations to perform.
#' It is only used when computing cell-specific differential expression between
#' groups specified in argument \code{data}.  
#' @param verbose logical that indicates if verbose messages should be shown.
#' 
#' @return Returns an NMF model object, with the cell specific signatures stored in the 
#' \code{\link{basis}} matrix, and the following miscellaneous slots that store details about the 
#' csSAM fit (accessible via \code{$}):
#' \item{deconv}{A list object containing a fit (cell-type specific
#' expression) for each group. Each element in the list is an object returned by
#' csFit.}
#' \item{fdr.csSAM}{ A list output of the fdrCsSAM function.}
#' \item{sigGene.csSAM}{A list of significant genes.}
#' 
#' @author 
#' Original function: Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
#' 
#' Adaptation for \pkg{CellMix}: Renaud Gaujoux
#' @seealso
#' \code{\link{csfit}},\code{\link{csSAM}},\code{\link{fdrCsSAM}}
#' @cite Shen-Orr2010
#' @import csSAM
.csSAM <- function(Y, x, data=NULL
					, nperms = 200, alternative = c('two.sided', 'greater', 'less')
					, standardize=TRUE, medianCenter=TRUE
					, logRm =FALSE, logBase = 2, nonNeg=TRUE
					, verbose = lverbose()){ 
		#function(G,cc,y,nperms = 200,alternative = 'two.sided'
#				,standardize=TRUE,medianCenter=TRUE, logRm =FALSE,logBase = 2,nonNeg=TRUE
#				,fileName='csSAMout.pdf') {
	
	## map arguments to those of csSamWrapper
	if( is.matrix(x) ){ # convert into an NMF model if necessary
		x <- nmfModel(Y, H=x, force.dim=FALSE)
	}
	
	# check and remove absent cell types
	if( any(nc <- apply(coef(x), 1L, function(x) all(x==0))) ){
		cn <- basisnames(x)[nc]
		wnote("Removed absent cell types from the estimation: "
				, str_out(cn), ' [', length(cn), ']')
		# subset on basis
		x <- x[!nc]
	}
	#
	
	if( isExpressionSet(Y) ) Y <- exprs(Y)
	# data contains the sample group definition
	if( is.null(data) ){
		data <- factor(rep(1, ncol(Y)))
	}else{
		if( !is.factor(data) ) data <- factor(data, levels=unique(data))
		# remove absent levels
		data <- droplevels(data)
	}
	# transpose/copy for csSAM
	G <- t(Y); cc <- t(coef(x)); y <- as.integer(data)
	##
	
	# remove NA-labeled samples
	if( length(i_rm <- which(is.na(y))) ){
		wnote('Dropped NA-labeled samples: ', str_out(rownames(G)[i_rm]), ' [', length(i_rm), ']')
		y <- y[-i_rm]
		cc <- cc[-i_rm, , drop=FALSE]
		G <- G[-i_rm, , drop=FALSE]
	}
	
	deconv <- list()
	
#	numset = length(unique(y))
	n <- summary(data, maxsum=Inf)
	numset <- nlevels(data)
	
	vmessage("Groups: ", str_out(n, Inf, use.names=TRUE, sep=' | '))
	# check parameters
	if( numset > 2L ){
		stop("csSAM - Cannot handle more than 2 groups of samples [", numset, "]")
	}
	if( nrow(G) != nrow(cc) ){
		stop("csSAM - Incompatible dimensions: number of cell proportions [", nrow(cc),"]"
			, " should match the number of samples [", nrow(G),"]")
	}
	#
	
	numgene = ncol(G)
	numcell = ncol(cc)
	cellNames = colnames(cc)
	
#	n <- vector(mode = "logical", length  = numset)
#	for (i in 1:numset) n[i] =sum(y==i)
	
#	for (curset in 1:numset) {
#		deconv[[curset]]= csfit(cc[y==curset,], G[y==curset,],logRm,logBase)
#	}
	vmessage("Fitting cell-specific linear model ... ", appendLF=FALSE)
	sets <- split(1:nrow(G), y)
	deconv <- lapply(sets, function(i){
				csfit(cc[i,,drop=FALSE], G[i,,drop=FALSE], logRm = logRm, logBase = logBase)
			})
	vmessage('OK')

	# early exit if only one group
	if( numset == 1L ){
		d <- dimnames(x)
		basis(x) <- t(deconv[[1L]]$ghat) 
		dimnames(x) <- d
		x$residuals <- t(deconv[[1L]]$residuals)
		x$se <- t(deconv[[1L]]$se)
		x$rawfit <- deconv[[1L]]
		return(x)
	}

	#rhat <- array(dim = c(numcell,numgene))
	vmessage("Computing csSAM model statistics ... ", appendLF=FALSE)
	rhat <- csSAM(deconv[[1]]$ghat, deconv[[1]]$se,
			n[1], deconv[[2]]$ghat, deconv[[2]]$se, n[2],
			standardize, medianCenter, nonNeg)
	vmessage("OK")
#	tt.sam <- runSAM(G, y)
	alternative <- match.arg(alternative)
	vmessage("Computing fdr using ", nperms, " permutations ... ", appendLF=FALSE)
	elapsed <- system.time({
	fdr.csSAM <- fdrCsSAM(G,cc,y,n,numcell,numgene, rhat
			, nperms = nperms, alternative = alternative
			, standardize = standardize, medianCenter = medianCenter
			, logRm = logRm, logBase = logBase, nonNeg = nonNeg)
	fdr.csSAM$alternative <- alternative
	})
	vmessage("OK")
	lmessage(2L, "Timing:")
	if( verbose >= 2L ) show(elapsed)
	
#	fdr.sam <- fdrSAM(G, y, nperms=nperms, tt.sam, alternative)
	vmessage("Finding signature genes ... ", appendLF=FALSE)
	sigGene <- findSigGene(G, cc, y, rhat, fdr.csSAM)
	vmessage('OK')
	
#	plotCsSAM(fdr.csSAM, fdr.sam,alternative,cellNames,numcell, file = fileName)
#	return(list(deconv = deconv, fdr.csSAM = fdr.csSAM, fdr.SAM = fdr.sam,sigGene.csSAM = sigGene,fileName = fileName))
	
	# wrap the result into the NMF model
	# remove permuation results if wrapping (to avoid memory issue)
	if( .cellmix.getOption('wrap') ) fdr.csSAM$rhatperm <- NULL
	x$rawfit <- list(deconv = deconv, fdr.csSAM = fdr.csSAM, sigGene.csSAM = sigGene
#					, fdr.SAM = fdr.sam
					)
	b <- t(rhat)
	rownames(b) <- rownames(Y)
	colnames(b) <- basisnames(x)
	basis(x) <- b
	x$call <- match.call()
	x$nperms <- nperms
	x$data <- data
	
	# return updated NMF model
	x
}



# overload csSAM function: use fast C++ code
# Will eventually be included in csSAM
#.findSigGene <- function (G, cc, y, rhat, csSAMData) 
#{
#	.Call('findSigGenes', rhat, csSAMData$cutp.g, csSAMData$fdr.g, PACKAGE='CellMix')
#}

# plain R version of findSigGene to test improvement in speed
.findSigGene_R <-
		function(G,cc,y,rhat,csSAMData) {
	numgene=ncol(G)
	numcell = ncol(cc)
	thresholdVec = csSAMData$fdr.g
	thresholdLen = length(thresholdVec[numcell,])
	sigGene <- array(dim = c(numcell, numgene))
	sigGene[,] = 1
	
	for (curThresh in 1:thresholdLen) {
		for (curcell in 1:numcell) {
			for (curgene in 1:numgene) {
				if(abs(rhat[curcell,curgene]) >= abs(csSAMData$cutp.g[curcell,curThresh])) {
					sigGene[curcell,curgene] = csSAMData$fdr.g[curcell,curThresh]
				}
			}
		}
	}
	
	return (sigGene)
}

#' Partial Gene Expression Deconvolution with csSAM
#' 
#' Estimates cell/tissue proportions given a known set of cell/tissue-specific 
#' expression signatures, using standard least-squares, as implemented 
#' by the package \code{\link{csSAM}}.
#'
#' All regressions are fitted using the function \code{\link{lsfit}}.
#'  
#' @inheritParams .csSAM
#' @aliases csSAM-ged
#' @cite Shen-Orr2010
#' @examples 
#' 
#' # random global expression
#' x <- rmix(3, 100, 20)
#' basisnames(x) <- paste('Cell', 1:nbasis(x))
#' # extract true proportions
#' p <- coef(x)
#' 
#' # deconvolve using csSAM
#' res <- ged(x, p, 'csSAM')
#' head(basis(res))
#' # proportions are not updated
#' identical(coef(res), p)
#' \dontshow{ 
#' 	stopifnot(identical(coef(res), p))
#'	stopifnot( nmf.equal(res, ged(x, p, 'csSAM')) ) 
#' }
#' 
#' # estimate cell-specific differential expression between 2 groups
#' gr <- gl(2, 10)
#' res <- ged(x, p, 'csSAM', data = gr, nperms=20, verbose=TRUE)
#' head(basis(res))
#' # plot FDRs
#' csplot(res)
#' # extract fdr for top differentially expressed gene in each cell type
#' t <- csTopTable(res)
#' str(t)
#' 
gedAlgorithm.csSAM <- setGEDMethod(key='csSAM'
		, description = "Estimates cell/tissue specific signatures from known proportions using SAM"
		, algorithm = .csSAM
		, reqBasis = FALSE
		, reqCoef= TRUE
		, reqMarker = FALSE
		, maxIter = 1L
		, cite = "Shen-Orr2010" 
)

#' The S3 method \code{csTopTable} for csSAM fits returns, for each feature, the false discovery 
#' rates of differential expression between groups of samples within each cell type, 
#' as computed by \code{\link[=csSAM]{fdrCsSAM}} when running csSAM.
#' These are returned as a list, whith one element per cell type.
#' 
#' @seealso \code{\link[csSAM]{fdrCsSAM}}, \code{\link{csTopTable}}
#' 
#' @inheritParams csTopTable.matrix
#' 
#' @S3method csTopTable ged_csSAM
#' @rdname gedAlgorithm.csSAM
csTopTable.ged_csSAM <- function(x, ...){
	
	# extract fitted model
	if( is.list(x) && length(x) && isNMFfit(x[[1L]]) ){
		x <- fit(x[[1L]])
		csSAMdata <- x$rawfit
		if( is.null(csSAMdata$sigGene.csSAM) ){
			stop("Cannot compute top table: csSAM fit does not contain gene significance data.\n  Was csSAM run with a group variable?")
		}
		res <- t(csSAMdata$sigGene.csSAM)
		rownames(res) <- rownames(x)
		colnames(res) <- basisnames(x)
	}else{
		if( is.null(x$sigGene.csSAM) ){
			stop("Cannot compute top table: csSAM fit does not contain gene significance data.\n  Was csSAM run with a group variable?")
		}
		res <- t(x$sigGene.csSAM)
	}
	
	csTopTable(res, ...)
}

#' The S3 method \code{csplot} for csSAM fits plots cell-specific fdr cumulative distributions. 
#' 
#' @param types index or names of the type to plot.
#' They need to be found in the fit data. 
#' @inheritParams graphics::plot.default 
#' 
#' @S3method csplot ged_csSAM
#' @rdname gedAlgorithm.csSAM
csplot.ged_csSAM <- function(x, types=NULL, xlab='# called', ylab='FDR', ylim=c(0,1), ...){
	
	# extract fitted model
	if( is.list(x) && length(x) && isNMFfit(x[[1L]]) ){
		model <- fit(x[[1L]])
		csSAMdata <- model$rawfit$fdr.csSAM
	}else{
		csSAMdata <- x
	}
	
	if( is.null(csSAMdata$rhat) ){
		warning("Nothing to plot: data does not contain any fdr data (did csSAM run with a group variable in argument `data`?)")
		return(invisible())
	}
	
	ncell <- nrow(csSAMdata$rhat)
	cellnames <- rownames(csSAMdata$rhat)
	
	if( is.null(cellnames) ){
		if( is.character(types) ){
			if( length(types) > ncell ){
				stop("Invalid number of types (", length(types), "):"
						, " should have at most the same number of cell type data in `x` (", ncell, ")")
			}
			cellnames <- types
		}else cellnames <- as.character(seq(ncell))
	}
	if( is.null(types) ) types <- seq(ncell)
	if( is.character(types) ){ # partial match the types
		types <- pmatch(types, cellnames)
		types <- types[!is.na(types)]
		if( !length(types) )
			stop("None of the types ", str_out(types), " matched [available: ", str_out(cellnames),"]")
	}
	
#	pdf(file = fileName, height = 11, width = 8)
	op <- par(mfrow = mfrow(ncell))
	on.exit( par(op) )
	
#	plot(SAMdata$ncall.sam, SAMdata$fdr.sam, xlab = "## called", 
#			ylab = "FDR", type = "l", log = "x", ylim = c(0, 1))
#	title(paste("SAM", alternative))
	
	lapply(types, function(i){
				j <- which(csSAMdata$ncall.g[i, ]<=0)
				plot(csSAMdata$ncall.g[i, -j], csSAMdata$fdr.g[i, -j]
						, xlab = xlab, ylab = ylab
						, type = "l"
						, log = "x", ylim = ylim)
				title(str_c(cellnames[i], ' - ', csSAMdata$alternative))
			})
	invisible(csSAMdata)
#	dev.off()
	
}
