# Test of algorithms for deconvolution
# 
# Author: Renaud Gaujoux
# Creation: 09 Feb 2012
###############################################################################


test.quantile <- function(n=1, quant=TRUE, rng){
	
	# load data
	curated <- readRDS('~/projects/hetero/paper-meegid/data/curated.RDS')	
	X <- exprs(curated$data$eset)
	# update markers
	m <- curated$data$markers
	m <- MarkerList(setNames(m@.Data, names(m)))
	# normalise
	Xn <- X
	if( quant ){
		Xn <- normalize.quantiles(X)
		rownames(Xn) <- rownames(X)
	}
	
	# perform semi-supervised
	set.seed(123456)
	res <- ged(Xn, m[,1:7], method='ssBrunet', nrun=n, .opt='v', seed='rprop')
	
	# plot results
	op <- par(mfrow=c(1,2))
	on.exit(par(op))
	profplot(curated$data$model, res, scale=TRUE)
	basismarkermap(m, res)
	res
	
}

test.alm <- function(n=1, method='flex', maxIter=100, nm=1, quant=FALSE, ...){
	
	# load data
	curated <- readRDS('~/projects/hetero/paper-meegid/data/curated.RDS')	
	X <- exprs(curated$data$eset)
	# update markers
	m <- curated$data$markers
	m <- MarkerList(setNames(m@.Data, names(m)))
	
	if( quant ){
		Xn <- normalize.quantiles(X)
		rownames(Xn) <- rownames(X)
		X <- Xn
	}
	
	# perform semi-supervised
	set.seed(12345)
	res <- nmf(X, 4, markers=m[, 1:nm], method=method, nrun=n, .opt='tv2', maxIter=maxIter, ...)
	
	# plot results
	op <- par(mfrow=c(1,2))
	on.exit(par(op))
	profplot(curated$data$model, res, scale=TRUE)
	basismarkermap(m, res)
	res
	
}
