# GED with the deconf method from Repsilber et al. (2010)
#
# Author: Renaud Gaujoux
###############################################################################

#' @include algorithm-common.R
NULL

# add extra post-installation action
setPackageExtra('install.packages', 'deconf', pkgs='deconf', repos='http://web.cbio.uct.ac.za/~renaud/CRAN')

###############################################################################
# Deconfounding approach from Repsilber et al (2010).
###############################################################################

# fast implementation of deconf
deconfounding2 <- function(target, x, maxIter = 1000L, data=NULL, error.threshold = 0){
	
	# show which version runs
	vmessage("Algorithm version 'fast'")
	
	# build call to .gedlab
	ca <- match.call()
	ca[['cscale']] <- TRUE
	ca[['sscale']] <- TRUE
	ca[['markers']] <- 'post'
	ca[[1L]] <- as.name('.gedlab')
	
	# run
	e <- parent.frame()
	eval(ca, e)
}

# Registration of NMF method 'deconf' [fast version]
nmfAlgorithm.deconf <- setNMFMethod('deconf', deconfounding2
		, defaults = list(seed='rprop')
		# define the objective function (infinite norm for matrices = larger eigen value)
		, objective = function(x, target, ...){ norm(fitted(x)-target, '2') }
		, overwrite = TRUE)

# Registration of NMF method '.deconf' [original version of 'deconf']
nmfAlgorithm._deconf <- setNMFMethod('.deconf', 'deconf'
		, defaults = list(seed = 'none')
		, algorithm=function(target, x, maxIter=1000L, data=NULL, error.threshold=0, ...){
			
			# show which version runs
			vmessage("Algorithm version 'original'")
			
			# use original implementation from the deconf package
			if( require.quiet('deconf', character.only=TRUE) ){
				deconfounding <- get('deconfounding', asNamespace('deconf'))
			}else{
				# try library where NMF in dev mode
				if( !isDevNamespace() ){
					stop("The package 'deconf' is required to run the original implementation of this algorithm.")
				} else{
					lib <- dirname(path.package('CellMix', quiet=TRUE))
					if( !require.quiet('deconf', character.only=TRUE, lib=file.path(lib, 'lib')) )
						stop("The package 'deconf' is required to run the original implementation of this algorithm.")
				}
			}
			
			# warn that the initial seed is not used
			if( isNumber(x) ) r <- x
			else if( is.nmf(x) ){
				r <- nbasis(x)
				if( !is.empty.nmf(x) )
					warning("Discarding NMF seed: running original 'deconf' method with ", r, " components [seed is generated in deconf::deconfounding].", immediate.=TRUE)
			}else{
				stop('Invalid seed/rank specification: must be either a number of an NMF model [', class(x), ']')
			}
			# apply the deconfounding algorithm
			if( nmf.getOption('verbose') > 1 )
				res <- deconfounding(target, n.cell.types=r, n.iterations = maxIter, error.threshold = error.threshold)
			else
				capture.output(res <- deconfounding(target, n.cell.types=r, n.iterations = maxIter, error.threshold = error.threshold))
			
			# wrap up the result
			basis(x) <- res$S$Matrix
			coef(x) <- res$C$Matrix
			residuals(x) <- res$error
			niter(x) <- res$nsim
			
			# assign components to cell-types using markers
			if( !is.null(data) ){
				vmessage("Assigning estimated components to cell-types using markers ... ", appendLF=FALSE)
				x <- match.nmf(x, data)
				vmessage('OK')
			}
			
			x
		}
		, overwrite = TRUE
)

#' Complete Gene Expression Deconvolution: Method deconf
#' 
#' The \code{\link{ged}} method \sQuote{deconf} uses an alternate least-squares 
#' algorithm to estimate both cell proportions and cell-specific signatures 
#' from global expression data, as proposed by \cite{Repsilber2010}.
#' 
#' This method fits an NMF model to the data in a completely \emph{unsupervised} manner. 
#' If marker genes are provided, they are used \strong{a posteriori} to assign each estimated 
#' component, i.e. each cell-specific signature, to the cell-type with the greatest proportions 
#' of consistent markers.
#' 
#' @section Fast built-in implementation:
#' The method \sQuote{deconf} is implemented as an NMF algorithm, which is registered 
#' under the same names in the \pkg{NMF} package's algorithm registry.
#' 
#' It uses an improved implementation, based on the fast combinatorial nonnegative least-squares 
#' algorithm from \cite{VanBenthem2004,KimH2007}, as provided by the 
#' function \code{\link{fcnnls}} in the \pkg{NMF} package.  
#' This enables to achieve great performance speed-up, being really -- way -- 
#' much faster than the original implementation.
#'
#' @aliases deconf-ged
#' @inheritParams .gedlab
#' @inheritParams deconf::deconfounding
gedAlgorithm.deconf <- setGEDMethod(key='deconf'
		, description = "Alternate least-square NMF method, using heuristic constraints [fast]"
		, algorithm = 'deconf'
		, reqBasis = FALSE
		, reqCoef= FALSE
		, reqMarker = FALSE
		, maxIter = 1000L
		, cite = "Repsilber2010"
)

#' @section Original implementation:
#' The \pkg{CellMix} also includes a way to run the original version from the \pkg{deconf} package, 
#' registered under the name \sQuote{.deconf} in the algorithm registry.
#' 
#' This version requires the \pkg{deconf} package, which was released as 
#' supplementary data only to support the paper from \cite{Repsilber2010}, i.e. it is not 
#' available from CRAN or Bioconductor.
#' However, we made it available from the \pkg{CellMix} CRAN-like support repository:
#' 
#' \url{http://web.cbio.uct.ac.za/~renaud/CRAN}
#' 
#' The easiest way to install it is to run:
#' 
#' \code{install.extras('CellMix', 'deconf')}
#' 
#' @rdname gedAlgorithm.deconf
#' @aliases .deconf-ged
gedAlgorithm._deconf <- setGEDMethod(key='.deconf'
		, description = "Alternate least-square NMF method, using heuristic constraints [original]"
		, algorithm = '.deconf'
		, reqBasis = FALSE
		, reqCoef= FALSE
		, reqMarker = FALSE
		, maxIter = 1000L
		, cite = "Repsilber2010"
)

#' Gene Expression Deconvolution: Flexible Experimental Algorithm 
#' 
#' \code{.gedlab} is an internal function of \pkg{CellMix} that implements a 
#' flexible workhorse function, whose arguments can be tuned to generate 
#' different kind of ged methods.
#' Its purpose is mainly to experiment and benchmark combinations of techniques, 
#' in order to design "better" deconvolution algorithm.
#' For example, it is called by the \sQuote{deconf} algorithm (fast version).  
#' 
#' @param v target matrix
#' @param x factorisation rank, i.e. the number of cell types to extract, or a 
#' complete initial NMF model.
#' @param data optional marker list used to either \emph{a posteriori} assign 
#' estimated signatures, or \emph{a priori} enforce marker block patterns 
#' \strong{before} each iteration.
#' @param maxIter maximum number of iterations
#' @param markers indicates what the markers are used for:
#' \describe{
#' \item{sQuote{prior}}{uses \code{\link[=DSA-ged]{DSA}} proportion estimation 
#' method from \cite{Zhong2013} to compute sensible initial proportions from 
#' average marker expression profiles in the mixed sample data.}
#' \item{sQuote{semi}}{enforces marker block patterns after each iteration.}
#' \item{sQuote{post}}{\emph{a posteriori} assigns estimated signatures;}
#' }
#' @param cscale logical that indicates if the estimated coefficient matrix, 
#' i.e. the proportions, should be scaled to force them to sum up to one, after 
#' each of its update.
#' @param sscale logical that indicates if the estimated signature matrix, 
#' i.e. the cell-specific expression profiles, should be scaled to force 
#' them to have a unit euclidean norm, after each of its update.
#' If \code{sscale=2}, then an inverse scaling is also applied to the 
#' coefficient matrix, so that the overall deviance of the model is not changed by 
#' this transformation.
#' @param ratio a numeric value that indicates the maximum ratio allowed between 
#' marker expression values in their own cell type and in other cell types.
#' It is meaningful only if greater than 1, but no errors is thrown if it is not.  
#' E.g. using \code{ratio=2} means that marker genes will have their expression 
#' values on signatures other than their own forced to a value at least twice lower 
#' than on their own cell types.
#' If \code{NULL}, then marker expression on cell types other than their own is 
#' forced to zero.
#' @param ... extra parameters currently not used
#' 
#' @keywords internal
#' @export
.gedlab <- function(target, x, data=NULL, maxIter=1000L
					, markers=c('prior+semi', 'prior', 'semi', 'post')
					, cscale=FALSE, sscale=FALSE
					, ratio=2, ...){

	verbose <- nmf.getOption('verbose')
	
	# enforce marker patterns if provided
	doAssign <- doEnforce <- doSeed <- FALSE
	M <- NULL
	if( !is.null(data) ){
		M <- matchIndex(data, target)
		if( !all(sapply(M, is.integer)) ){
			stop('Unexpected error - invalid marker list: matched marker IDs must be integer indexes.')
		}
		markers <- match.arg(markers)
		vmessage("Using ", nmark(M), "/", nmark(data), " markers [", markers, ']')
		lmessage(2, paste0(" ", capture.output(print(nmark(M, TRUE))), collapse="\n"))
		
		# define how to use markers 
		doAssign <- grepl('post', markers)
		doEnforce <- grepl('semi', markers)
		doSeed <- grepl('prior', markers)
	}
	#

	v <- target
	# static variables
	error.W <- Inf
	error.H <- Inf
	error <- Inf
	min.error <- Inf
	
	# show scaling strategy
	lmessage(2, "Scaling Strategy - coef: "
				, if(cscale) 'yes' else 'no'
				, ' | signatures: '
				, if(sscale) 'yes' else 'no'
				, if(sscale>1) ' (symetric)' )
	
	# initial column normalization
	I.noised <- if( sscale ) .apply.constraints.S(v) else v
	
	if( doSeed ){# seed from markers: DSA + fcnnls
		doAssign <- FALSE
		if( verbose ) message("Seeding proportions using DSA method ... ", appendLF=FALSE)
		coef(x) <- .DSAproportions(v, M)
		if( verbose ) message("OK")
	}else{
		r <- NULL
		if( is.empty.nmf(x) ){ # use rank from empty model
			message('# NOTE: Initialise empty NMF model internally (rprop)')
			r <- nbasis(x)
		}else if( isNumber(x) ){  
			# initialise an NMF model if only the rank is provided
			message('# NOTE: Initialise NMF model internally (rprop)')
			r <- x
		}else if( !is.nmf(x) ){
			stop("Invalid model/seed specification: expected an NMF model or a number [", class(x), ']')
		}
		# seed with rprop if necessary
		if( !is.null(r) ) x <- seed(r, I.noised, method='rprop')
	}
	
	# seed basis signature enforcing markers 
	if( doEnforce ){
		if( verbose ) message("Computing initial cell-specific signatures using qprog/fcnnls method ... ", appendLF=FALSE)
		basis(x) <- csfit_markers(v, coef(x), M, ratio=ratio)
		if( verbose ) message("OK")
	}

	# alternating: either take S as fix and calculate C or take C
	# as fix and calculate S.
	# start: S is considered as fix
	
	# apply constraints to basis and coef
	if( sscale ) basis(x) <- .apply.constraints.S(basis(x))
	if( cscale ) coef(x) <- scoef(x)
	
	nrestart <- 0L
	i <- 0L
	
	if( verbose ) pgr <- iterCount(maxIter)
	stop.signal <- FALSE
	while( i < maxIter){
		i <- i + 1L
		if( verbose ) pgr(i)
		# force marker expression patterns on signatures
		if( doEnforce ){
			.basis(x) <- neq.constraints.inplace(basis(x), M, ratio=ratio, copy=TRUE)
		}
		#
	
		# nnls for H
		if( i>1 || !doSeed ){
			nlsfit <- .lsqnonneg(basis(x), I.noised, check=FALSE)
			if (	nlsfit$resnorm >= error || nlsfit$resnorm == error.H  ) #no changes, we stuck in local minima
			{
				if( verbose > 1 ) cat(" - Found minimum on H [", nlsfit$resnorm, "].\n", sep='')
				stop.signal <- TRUE
				break;
			}
			error <- error.H <- nlsfit$resnorm
			if ( min.error > error ) min.error <- error
			
			# store in x
			.coef(x) <- nlsfit$x
		
			if ( any(zr <- rowSums(nlsfit$x)==0) ){
				if( verbose > 1 ) cat(" - ", sum(zr), " null row(s) in H => restarting... \n");
				nrestart <- nrestart + 1L;
				if ( nrestart >= 10L ){
					warning("Too many ALS restarts due [Computation stopped after the 9th restart]");
					stop.signal <- TRUE
					break;
				}
				
				# re-initialize random W
				# re-initialize base average error
				error.W <- Inf
				error.H <- Inf
				error <- Inf
				min.error <- Inf
				basis(x) <- .apply.constraints.S( rmatrix(basis(x)) )
				next;
			}
			
			# rescale as proportions
			if( cscale ) .coef(x) <- scoef(x)
		}
		
		# nnls for W
		nlsfit <- .lsqnonneg(t(coef(x)), t(I.noised), check=FALSE)
	
		if (	nlsfit$resnorm >= error || nlsfit$resnorm == error.W ) #no changes, we stuck in local minima
		{
			if( verbose > 1 ) cat(" - Found minimum for W [", nlsfit$resnorm, "].\n", sep='')
			stop.signal <- TRUE
			break;
		}
		error <- error.W <- nlsfit$resnorm
		if ( min.error > error ) min.error <- error
		
		# store in x
		w <- t(nlsfit$x)
		.basis(x) <- w
		
		# rescale signatures
		if( sscale ){
			a <- colMeans(basis(x), na.rm=TRUE)
			.basis(x) <- .apply.constraints.S(basis(x))
			if( sscale > 1) # apply inverse transformation to H
				.coef(x) <- sweep(.coef(x), 1L, a, '*')
		}
		
		# track error (after transformation)
		x <- trackError(x, norm(fitted(x)-I.noised, '2'), niter=i)
		
		#show the error
		if( verbose > 2 ){
			lmessage(3, " - error: H= ", error.H, " | W= ", error.W, " - min.error= ", min.error)
		}
	}
	
	if( verbose ){
		ended <- if( stop.signal ) 'converged' else 'stopped'
		if( verbose > 1 ) cat("DONE")
		cat(" (", ended, " at ",i ,'/', maxIter," iterations)\n", sep='')
	}
	
	# set number of iterations performed
	niter(x) <- i
	
	# assign signatures using markers
	if( doAssign ){
		vmessage("Assigning estimated components to cell-types using markers ... ", appendLF=FALSE)
		x <- match.nmf(x, data)
		vmessage('OK')
	}
	#
	
	# put dimnames back in if necessary
	if( doAssign || doSeed || doEnforce ){
		if( is.null(basisnames(x)) ) basisnames(x) <- names(M)
	}
	
	x
}


nmfAlgorithm.gedlab <- setNMFMethod('.gedlab', .gedlab
		, defaults = list(seed='rprop')
		# define the objective function (infinite norm for matrices = larger eigen value)
		, objective = function(x, target, ...){ norm(fitted(x)-target, '2') }
		, overwrite = TRUE)

gedAlgorithm.hybrid <- setGEDMethod(key='.hybrid'
		, description = "Experimental hybrid algorithm"
		, algorithm = '.gedlab'
		, reqBasis = FALSE
		, reqCoef= FALSE
		, reqMarker = FALSE
		, maxIter = 1000L
)


########################
## UTILS
########################

.apply.constraints.S <- function(x){
	sweep(x, 2L, colMeans(x, na.rm=TRUE), '/')
}

