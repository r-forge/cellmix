# Simulation of global expression from mixed cell/tissue type
# 
# Author: Renaud Gaujoux
# Creation: 19 Dec 2011
###############################################################################


#' Generating Random Cell Type Proportions
#' 
#' \code{rproportions} generates a random NMF model, whose mixture coefficients
#' for each sample (i.e. column of the coefficient matrix) fulfil the sum-up-to-one 
#' constraint of proportions.
#'  
#' Internally, it generates a full random NMF model using the function \code{\link[NMF]{rnmf}}, 
#' and either scale the columns of its coefficient matrix so that they sum-up to one 
#' if \code{alpha=NA}, or re-draw them from a dirichlet distribution using the  
#' function \code{\link[gtools]{rdirichlet}} if \code{alpha} is a single number or a 
#' numeric vector.
#' 
#' This function is suitable for seeding NMF computation, and is in fact registered
#' as the seeding method \sQuote{rprop}, which can be used in calls to 
#' \code{\link[NMF]{nmf}}.
#'  
#' @inheritParams NMF::rnmf
#' @param ... extra arguments passed to \code{\link[NMF]{rnmf}}.
#' @param alpha shape parameter(s) for the dirichlet distribution.
#' It must be a single numeric value, or a numeric vector of length the 
#' number of basis components in the NMF model -- which is generally specified 
#' in \var{object}.
#' If a single numeric, it is used as a shape parameters for all components.  
#' If \code{alpha=NA}, then the coefficients are only scaled to sum-up to one.
#' 
#' Note that because it appears after arguments \code{...}, its name must be 
#' fully specified, or it will be passed to \code{rnmf}, while its default value 
#' used within \code{rproportions}.
#'  
#' @export
#' @importFrom gtools rdirichlet
#' @examples 
#' 
#' # target data
#' x <- rmatrix(20,10)
#' 
#' # generate models
#' rproportions(3, x)
#' rproportions(3, 30, 15)
#' 
#' # scaling only
#' rproportions(3, 20, 10, alpha=NA)
#' # full shape specification
#' rproportions(3, x, alpha=c(1,2,3))
#' 
#' # error if alpha has wrong length
#' try( rproportions(3, 20, 10, alpha=c(1,2,3,4)) )
#' 
rproportions <- function(x, target, ..., alpha=1){
	
	# create random NMF model
	x <- rnmf(x, target, ...)
	
	# create random mixture proportions
	if( !is_NA(alpha) ){ # from dirichlet distribution
		if( isNumber(alpha) ) alpha <- rep(alpha, nbasis(x))
		coef(x) <- t(rdirichlet(ncol(x), alpha=alpha))
	}else{
		# impose sum up to one constraint
		coef(x) <- scoef(x)
	}
	
	x 
}

# register the new seeding method
setNMFSeed('rprop', rproportions, overwrite = TRUE)


#' Generating Random Global Mixed Gene Expression Data 
#' 
#' The function \code{rmix} generates an \code{\linkS4class{ExpressionMix}} object,
#' composed of a given number of underlying cell types.
#' The amount of noise added to both the cell-specific signatures and the global expression 
#' values is customisable.  
#'
#' @param x number of true underlying cell types or a matrix containing the signatures themselves, 
#' i.e. cell-specific expression values for each feature. 
#' For convenience, it may also specify the markers to enforce on the signatures, as a vector or list 
#' of length > 1, in which case argument \code{markers} must be missing.  
#' @param n number of features, i.e. genes.
#' The argument is required if \code{x} specifies the number of signatures.
#' If \code{x} is provided as a matrix, then \code{n} is used to subset it (\code{x[n, ]}) 
#' before simulating the global expression data.
#' @param p number of samples
#' @param markers specification of the number of markers to enforce on each cell type signature.
#' This should be a value supported by \code{\link{enforceMarkers}}.
#' Markers enforcement may be disabled with \code{markers=NA}.
#' @param min minimum cell-specific expression value before adding noise and marker differential.
#' @param max maximum cell-specific expression value before adding noise and marker differential.
#' @param mfold fold change expected on cell-specific expression for marker genes
#' @param alpha parameter for the dirichlet distribution from which are drawn the mixture 
#' proportions, using \code{\link[gtools]{rdirichlet}}.
#' @param snoise parameters for the normal noise added to each true underlying signatures as
#' \eqn{x + N(\mu, \sigma)}.
#' @param gnoise parameters for the normal noise with inverse gamma variance added to each 
#' feature global expression profile as \eqn{e_{ij} + N(0, 1/\gamma_i)}.
#' @param ... extra arguments currently not used.
#' 
#' @return an \code{\linkS4class{ExpressionMix}} object, that contains the true underlying 
#' signatures and proportions stored as an NMF model.
#' 
#' @importFrom gtools rdirichlet
#' @export
#' @examples 
#' 
#' # 3 cell types, 100 features, 20 samples
#' rmix(3, 100, 20)
#' 
#' # from known signature matrix
#' s <- rmatrix(100, 5)
#' x <- rmix(s, p=20)
#' dim(x)
#' if( !isCRAN_timing() ){
#' aheatmap(x)
#' }
#' 
#' # markers are enforced on each true signature
#' x <- rmix(4, 50, 20, markers=6)
#' if( !isCRAN_timing() ){
#' basismap(x, Rowv=NA)
#' }
#' # or also
#' x <- rmix(1:4, 50, 20)
#' if( !isCRAN_timing() ){
#' basismap(x, Rowv=NA)
#' }
#' 
rmix <- function(x, n=100, p=20, markers=ceiling(nrow(x)/20), min=0, max=20, mfold=2, alpha=1
					, snoise=list(mean=0, sd=0.05), gnoise=list(shape=5, scale=1), ...){
	
	# extract expression data if necessary
	if( isExpressionSet(x) ) x <- exprs(x)
	
	x <- if( x_is_mat <- is.matrix(x) ){
		if( !missing(n) ){
			# use first n rows if single number
			if( isNumber(n) ) n <- 1:n
			x <- x[n, ]
		}
		x
	}else{
		# use x as markers
		if( length(x) > 1L && missing(markers) ){
			markers <- x
			x <- length(x)
		}
		if( isNumber(x) ){
			if( !isNumber(n) ) 
				stop("Invalid argument `n`: expecting single numeric value (number of target rows) [", class(n), ']')
			rmatrix(n, x, min=min, max=max)
		} else stop("Invalid argument `x`: must be a matrix or a single numeric value. [", class(x) , "]")
	}


	# create random mixture proportions
	if( isNumber(alpha) ) alpha <- rep(alpha, ncol(x))
	prop <- t(rdirichlet(p, alpha=alpha))

	# build true underlying NMF model
	model <- nmfModel(x, prop)
	
	if( isTRUE(markers) ) markers <- ceiling(nrow(model)/20)
	# define cell type names
	CLnames <- names(markers) %||% paste0('CL_', 1:nbasis(model))
	
	# add markers if necessary
	if( length(markers) && !is_NA(markers) ){
		basisnames(model) <- CLnames
		model0 <- model
		model <- enforceMarkers(model, markers, ratio=mfold, attach=TRUE)
		# update signature matrix
		x <- basis(model)
#		i <- marknames(getMarkers(x))
		#x[i,] <- x[i,] + rmatrix(x[i,], max=max(x[i,])/mfold)
	}

	# add biological variability to the signatures
	xN <- abs(x + do.call('rmatrix', c(list(x, rnorm), snoise)))
	
	# combine signatures and proportions
	y <- xN %*% prop
	# add global noise (with inverse Gamma variance)
	lambda <- do.call('rgamma', c(list(nrow(y)), gnoise))
	y <- t(sapply(1:nrow(y), function(i){
		y[i,] + rnorm(y[i,], mean=0, sd=1/lambda[i])			
	}))
	
	# warp into an ExpressionMix object
	res <- ExpressionMix(abs(y), composition=model)
	# set dimnames
	if( !x_is_mat ){
		rownames(res) <- paste0('gene_', 1:nrow(res))
		basisnames(res) <- CLnames
	}
	colnames(res) <- paste0('sample_', 1:ncol(res))
		
	res
}

#' Generating Random Pure Cell Type Sample
#'
#' \code{rpure} generates random expression data that simulates
#' pure cell type samples, that share most gene expression profile pattern.
#' 
#' @param x number of cell types
#' @param n number of features
#' @param p number of samples per cell type
#' @inheritParams stats::runif 
#' 
#' @return an \code{\link{ExpressionMix}} object that contains
#' the expression data, the pure basis signatures, and a phenotype
#' variable \code{'CellType'}, which indicates the cell type of 
#' each sample.
#' 
#' @export
#' @examples 
#' 
#' x <- rpure(3)
#' aheatmap(x, annCol=TRUE)
#' 
rpure <- function(x, n=100, p=20, min=0, max=20){
	
	# generate correlated cell-specific profiles
	pure <- rmatrix(n, x, min=min, max=max)	
	exp <- do.call('cbind', lapply(1:ncol(pure), function(i) pure[,i] + rmatrix(length(pure[,i]), p, dist=rnorm)))
	ct <- list(CellType = gl(x, p, labels=paste0('CL-', 1:3)))
	# make a common part
	nc <- ceiling(n/2) 
	com <- runif(nc, min=min, max=max)
	exp[seq(n-nc+1, n),] <- com + rmatrix(nc, ncol(exp), dist=rnorm)
	
	# wrap everything in an ExpressionMix object
	ExpressionMix(exp, composition=nmfModel(W=pure)
				, phenoData=AnnotatedDataFrame(data.frame(CellType=ct)))
}

