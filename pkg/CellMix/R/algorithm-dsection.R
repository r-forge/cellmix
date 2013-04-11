# R implementation of DSection, from the MATLAB code available at:
# http://www.cs.tut.fi/~erkkila2/software/dsection/
# Accessed on 31 Oct 2011
#
# Reference:
# Probabilistic analysis of gene expression measurements from heterogeneous tissues
# Timo Erkkilä et al.
# Bioinformatics (2010) 26 (20): 2571-2577. doi: 10.1093/bioinformatics/btq406
# 
# Author: Renaud Gaujoux
# Creation: 31 Oct 2011
###############################################################################
#' @include registry-algorithms.R
NULL

# add extra post-installation action
setPackageExtraHandler('octave', function(pkgs, repos='http://web.cbio.uct.ac.za/~renaud/CRAN/extras/CellMix/'){
	
	# load RcppOctave
	if( !require(RcppOctave) ){
		warning("The package RcppOctave is required to perform Octave package installation: "
				, "skipped installation of ", str_out(pkgs, Inf))
		return()
	}
	
	# get package index
	index <- tempfile('octave_repo')
	on.exit( unlink(index) )
	download.file(repos, destfile=index)
	idx <- readLines(index)
	m <- str_match_all(idx, "(([-a-zA-Z_]+)-[0-9.]+\\.tar\\.gz)")
	m <- m[sapply(m, length)>0]
	m <- sapply(m, function(p) setNames(p[1,2], p[1,3]) )
	
	# match packages
	i <- match(pkgs, names(m))
	if( length(ipb <- which(is.na(i))) ){
		warning("Could not find Octave package(s) ", str_out(pkgs[ipb]), " in repository '", repos, "'")
	}
	i <- i[!is.na(i)]
	if( !length(i) ) return()
	found <- m[i]
	
	# get package
	dir <- tempfile('octave_packages')
	dir.create(dir, recursive=TRUE)
	on.exit( unlink(dir, recursive=TRUE), add=TRUE)
	pkgfile <- file.path(dir, found)
	remote <- str_c(repos, found)
	mapply(download.file, remote, pkgfile)
	
	cmd <- "sudo octave --silent --eval \"pkg install "
	sapply(pkgfile, function(p){
		message("Installing Octave package ", p)
		cmd <- str_c(cmd, "'", p, "'\"")
		system(cmd)
		#RcppOctave::.O$pkg('install', p)
	})
})

setPackageExtra('octave', 'DSection', pkgs=c('io', 'general', 'miscellaneous', 'struct', 'optim', 'statistics'))


#' DSection Gene Expression Deconvolution Method
#' 
#' The \emph{DSection} algorithm performs gene expression deconvolution when priors 
#' on proportions are available, using a Markov Chain Monte Carlo approach 
#' \cite{Erkkila2010}.
#' 
#' In \pkg{CellMix}, this method is registered with the key \code{'DSection'},
#' and is can be applied to gene expression data via the function \code{\link{ged}}. 
#'   
#' @details
#' This function uses the \pkg{RcppOctave} package to run the original Matlab code 
#' in \emph{Octave}.
#' The documentation was extracted from the Matlab source file, that can be 
#' found in the \pkg{CellMix} package "scripts/DSection" subdirectory.
#' 
#' The Matlab code requires the Octave packages \emph{statistics} and \emph{optim} to run properly.
#' These packages can be downloaded from Octave-forge: 
#' 
#' \url{http://sourceforge.net/projects/octave/files/Octave\%20Forge\%20Packages/Individual\%20Package\%20Releases/}
#' 
#' and installed in Octave with: 
#' 
#' \code{pkg install '<path/to/package/tar/gz/file>'}
#' 
#' or in R
#' 
#' \code{install.extras('CellMix', 'octave:DSection')}  
#' 
#' @param Y <I-by-J> matrix of measurements from heterogeneous tissues.
#' I is the number of probes/genes/etc., and J is the number of tissues.
#' @param p0 <T-by-J> matrix of prior predictions on cell type proportions.
#' T is the number of cell types, and columns in \code{p0} must be positive and 
#' add up to one.
#' @param groups <1-by-J> vector of treatment indices, so that 
#' unique(Treatment) = [1,2,...,C], where C is the number of treatments 
#' including control, i.e., "no treatment", if available.
#' @param W0 Prior prediction weight, i.e., degree of confidence, on \code{p0}.
#' Defines the peakedness of Dirichlet density around p0. NOTE: keep W0 >= T.
#' @param W_proposal Transition kernel weight, defines the peakedness of 
#' Dirichlet density around p*, the old value. The higher W_proposal is, 
#' the smaller the proposal steps around p* are. 
#' @param nBurnIn Amount of burn-in. NOTE: keep nBurnIn > 0.
#' @param nSamples Amount of sampling. NOTE: keep nSamples > 0.
#' @param samplep logical value, indicating whether to sample from the 
#' posterior for cell type proportions (\code{TRUE}) or not (\code{FALSE}).
#' SUGGESTED USE: sample from the posterior (samplep = 1).
#' @param summarize logical indicating whether only average values should be 
#' returned -- and computed.
#' @param verbose logical that indicates if verbose messages should be shown.
#' 
#' @return A list with the following elements:
#' 
#' \item{MCData}{ results from the MCMC estimation.}
#' \item{x_LS}{ Standard-least square estimate.}
#' \item{groups}{ factor defining the groups of samples, if any was provided.}
#' \item{call}{ the call to DSection.}
#' \item{parameters}{ a list of some of the parameters used in the estimation.}
#' \item{p0}{ initial prior on proportions.}
#'   
#' @family dsection
#' @source \url{http://www.cs.tut.fi/~erkkila2/software/dsection/DSection.m}
#' @author
#' Original Matlab code: Timo Erkkila 
#'  
#' Wrapper function: Renaud Gaujoux
#' @cite Erkkila2010
#' @references
#' \url{http://informatics.systemsbiology.net/DSection/}
#' \url{http://www.cs.tut.fi/~erkkila2/software/dsection/index.html}
#' 
#' @export
#' @examples 
#' 
#' #' # random global expression
#' x <- rmix(3, 50, 10)
#' dim(x)
#' 
#' # extract true proportions
#' p <- coef(x)
#' # add noise to proportions
#' p0 <- scoef(abs(p + rmatrix(p, dist=rnorm, sd=0.15)))
#' # check how noisy this makes the proportion estimates
#' profplot(p, p0)
#' 
#' \dontrun{
#' # Requires some extra Octave packages to be installed (it defines `gamfit`)
#' res <- DSection(exprs(x), p0, nSamples=100, nBurnIn=1L)
#' profplot(p, res$MCData$p)
#' }
#' 
DSection <- function(Y, p0, nSamples, groups=NULL, nBurnIn=4*nSamples
						, W0=max(nrow(p0), 10), W_proposal=100, samplep=TRUE, summarize=TRUE
						, verbose=FALSE){
	
	# load RcppOctave
	if( !require(RcppOctave) )
		stop("The package RcppOctave is required to run this function.")
	
	# source required m-files
	o_source(packagePath('scripts', 'DSection', 'DSection.m'))
	o_source(packagePath('scripts', 'DSection', 'MCMCSummary.m'))
	
	Treatment <- groups
	if( is.null(Treatment) ) Treatment <- rep(1,ncol(Y))
	# check arguments
	if( !is.factor(Treatment) ) Treatment <- factor(Treatment)
	
	# call Matlab function
	.CallOctave('pkg', 'load', 'optim', argout=0)
	res <- .CallOctave('DSection', Y, p0, W0, W_proposal, as.numeric(Treatment), nBurnIn, nSamples, samplep, summarize, verbose)
	# Y,p0,W0,W_proposal,Treatment,nBurnIn,nSamples,samplep,summarize,verbose
	
	# transpose results
	if( nlevels(Treatment) == 1L ){
		res$MCData$x <- t(res$MCData$x)
		res$x_LS <- t(res$x_LS)
		dn <- list(rownames(Y), rownames(p0))
	}else{
		x <- array(dim=c(nrow(Y), nrow(p0), nlevels(Treatment)))
		x_LS <- array(dim=c(nrow(Y), nrow(p0), nlevels(Treatment)))
		sapply(1:nlevels(Treatment), function(i){
					x[,,i] <<- t(res$MCData$x[,,i])
					x_LS[,,i] <<- t(res$x_LS[,,i])
				})
		dn <- list(rownames(Y), rownames(p0), levels(Treatment))
		res$MCData$x <- x
		res$x_LS <- x_LS
		res$groups <- Treatment
		res$MCData$groups <- Treatment
	}
	# set dimnames
	if( !all(vapply(dn, is.null, logical(1))) ){
		dimnames(res$MCData$x) <- dn
		dimnames(res$x_LS) <- dn
	}
	
	# add extra parameters
	res$call <- match.call()
	res$parameters <- list()
	res$parameters$nSamples <- nSamples
	res$parameters$nBurnIn <- nBurnIn
	res$p0 <- p0
	res$parameters$W0 <- W0
	res$parameters$W_proposal <- W_proposal
	res$parameters$Treatment <- Treatment
	res$parameters$samplep <- samplep
	res$parameters$summarize <- summarize
	
	# return result
	res
}

# Registration of NMF method 'DSection'
nmfAlgorithm.DSection <- setNMFMethod('DSection',
	function(target, seed, maxIter, p0, data=NULL, nBurnIn=4*maxIter
			, W0=max(100, nbasis(seed)), W_proposal=100, samplep=TRUE, ...){
		
		# call DSection with nSamples = maxIter
		verbose <- nmf.getOption('verbose')
		if( nmf.getOption('verbose') > 1 ){
			res <- DSection(target, coef(seed), nSamples = maxIter, groups=data
							, nBurnIn=nBurnIn, W0=W0, W_proposal=W_proposal, samplep=samplep
							, summarize=TRUE, verbose=verbose, ...)
		}else{
			capture.output(res <- DSection(target, coef(seed), nSamples = maxIter, groups = data
								, nBurnIn=nBurnIn, W0=W0, W_proposal=W_proposal, samplep=samplep 
								, summarize=TRUE, verbose=verbose, ...))
		}
		
		## WRAP UP
		# basis
		if( is.null(data) ){
			basis(seed) <- res$MCData$x
		}else{
			# store basis for first group
			basis(seed) <- res$MCData$x[,,1]
		}
		# store whole result in misc slot
		seed$rawfit <- res$MCData
		# coef
		coef(seed) <- res$MCData$p
		# number of iteration = number of recorded sampling
		niter(seed) <- maxIter
		# add miscellaneous results
		seed$lambda <- res$MCData$lambda
		if( !is.null(res$groups) )
			seed$groups <- res$groups
		if( !is.null(res$x_LS) )
			seed$x_LS <- res$x_LS
		if( !is.null(res$lambda_LS) )
			seed$lambda_LS <- res$lambda_LS
		# add extra parameters
		seed@parameters <- mergeList(seed@parameters, res$parameters)		
		##
		
		# TODO: compute statistics
		
		# return updated seed
		seed
	} 
	, objective = "euclidean"
	, overwrite = TRUE
)	

#' Partial Gene Expression Deconvolution with DSection
#' 
#' Estimates cell/tissue cell/tissue-specific expression signatures, 
#' given proportion priors using the MCMC approach from \cite{Erkkila2010}, 
#' implemented in Matlab and wrapped in \pkg{CellMix} by the function 
#' \code{\link{DSection}}.
#'
#' @inheritParams DSection
#' @param data variable (e.g., factor) that defines groups of samples.
#' Cell-specific signatures will be computed for each group.
#' @param maxIter number of sampling, not including burn-in sampling.
#' @param ... extra arguments passed to \code{\link{DSection}}.
#' @aliases DSection-ged
#' @examples 
#' 
#' # random global expression
#' x <- rmix(3, 50, 20)
#' # extract true proportions
#' p <- coef(x)
#' # add noise
#' p0 <- scoef(abs(p + rmatrix(p, dist=rnorm, sd=0.15)))
#' # check how noisy this makes the proportion estimates
#' profplot(p, p0)
#' 
#' # deconvolve using DSection
#' \dontrun{
#' res <- ged(x, p0, 'DSection', maxIter=10, seed=12345)
#' head(basis(res))
#' # proportions are updated
#' !identical(coef(res), p0)
#' # check how better they are
#' profplot(x, res)
#' 
#' \dontshow{ 
#' 	stopifnot( !identical(coef(res), p0) )
#'	stopifnot( nmf.equal(res, ged(x, p0, 'DSection', maxIter=10, seed=12345)) ) 
#' }
#' }
#' 
#'
gedAlgorithm.DSection <- setGEDMethod(key='DSection'
, description = "Estimates proportions from proportions priors using MCMC" 
, algorithm = 'DSection'
, reqBasis = FALSE
, reqCoef= TRUE
, reqMarker = FALSE
, maxIter = 500L
, cite = "Erkkila2010"
)


#' The S3 method \code{csTopTable} for DSection fits computes nominal p-values (i.e. unadjusted)
#' of differential expression between cell type or group of samples within each 
#' cell type, for deconvolution results from the \code{\link[=DSection-ged]{DSection}} algorithm.
#' 
#' @inheritParams csTopTable.matrix
#' @param coef specifies the reference cell type or group of samples, 
#' for which differential expression is computed. 
#'   
#' @S3method csTopTable ged_DSection
#' @rdname gedAlgorithm.DSection
csTopTable.ged_DSection <- function(x, coef=1L, decreasing=TRUE, ...){
	
	# extract fitted model
	if( is.list(x) && length(x) && isNMFfit(x[[1L]]) ){
		MCData <- x[[1L]]$rawfit
	}else{
		MCData <- x
	}
	
	dims <- dim(MCData$x)
	ncell <- ncol(MCData$x)
	ngenes <- nrow(MCData$x)
	cellnames <- colnames(MCData$x)
	lambda <- as.numeric(MCData$lambda)
	
	if( length(dims) == 2L ){
		nsamples <- ncol(MCData$p)
		res <- (MCData$x[, -coef, drop=FALSE] - MCData$x[, coef, drop=TRUE])
		res <- res / sqrt(1/MCData$lambda * 2/nsamples)
		# add dimnames
		colnames(res) <- cellnames
		rownames(res) <- rownames(MCData$x)
	}else{
		nsamples <- summary(MCData$groups, max.sum=Inf)
		res <- lapply(seq(length(nsamples))[-coef], 
				function(j){
					m <- (MCData$x[, , j, drop=FALSE] - MCData$x[, , coef, drop=FALSE])
					res <- m / sqrt(1/lambda * sum(1/nsamples[c(j, coef)])) 
					res[,,1]
				})
		if( length(res) == 1L ){
			res <- res[[1L]]
			# add dimnames
			colnames(res) <- cellnames
			rownames(res) <- rownames(MCData$x)
		}else{ # array
			a <- array(NA, dim=c(length(res), dim(res[[1L]])))
			lapply(seq_along(res), function(i){
						x <- res[[i]]
						a[i,,] <<- x
					})
			# add dimnames
			dimnames(a)[[3L]] <- cellnames
			dimnames(a)[[2L]] <- rownames(MCData$x)
			#
			res <- a
			dimnames(res)[[1L]] <- names(nsamples)[-coef]
		}
	}
	
	csTopTable(res, ..., decreasing=decreasing)
}

.test.DSection <- function(n=10, b=10, loaddata=FALSE, ...){
	
	# load data
	if( loaddata )
		load('~/projects/hetero/paper-meegid/data/curated.rda', .GlobalEnv)
	curated <- get('curated', .GlobalEnv)
	X <- exprs(curated$data$eset)
	X <- log2(X)
	
	r <- 4
	p0 <- matrix(0.25, r, ncol(X)) + rmatrix(r, ncol(X), rnorm, mean=0, sd=0.5)
	p0 <- coef(curated$data$model) + rmatrix(r, ncol(X), rnorm, mean=0, sd=0.1)
	op <- par(mfrow=c(1,2))
	on.exit(par(op))
	profplot(curated$data$model, p0)
	W0 <- 10
	W_proposal <- 100
	nBurnIn <- b
	nSamples <- n
	Treatment <- rep(1, ncol(X)) 
	samplep <- TRUE

	res <- DSection(X,p0,W0,W_proposal,Treatment,nBurnIn,nSamples,samplep,...)
	
	W <- t(res$MCData$x)
	H <- res$MCData$p
	rownames(W) <- rownames(X)
	res.model <- nmfModel(W, H)
	profplot(curated$data$model, res.model)
	res.model
}

