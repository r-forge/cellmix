# Unit test for GED algorithms
# 
# Author: Renaud Gaujoux
# Creation: 21 Jan 2012
###############################################################################



test.all <- function(){
	
	# retrieve all algorithms
	keys <- gedAlgorithm()
	
	# random data
	x <- rmatrix(100,20)
	model <- rnmf(3, x)
	marks <- gMarkerList(3, 5)
	maxIter <- 5 
	
	op <- cellmix.options(verbose=TRUE)
	on.exit( cellmix.options(op) )
	lapply(keys, function(k){
		
		if( k == 'DSection' && !require(RcppOctave) ) return(TRUE)
		
		# load algorithm
		meth <- GEDStrategy(k)
		
		# apply algorithm on relevant input data
		res <- 
		if( onlyRequired('Basis', meth) ) 
			checkTrue(isNMFfit(ged(x, basis(model), meth, maxIter=maxIter)), paste("Algorithm ", dQuote(k), " works with basis only."))
		else if( onlyRequired('Coef', meth) )
			checkTrue(isNMFfit(ged(x, coef(model), meth, maxIter=maxIter)), paste("Algorithm ", dQuote(k), " works with coef only."))
		else if( onlyRequired('Marker', meth) )
			checkTrue(isNMFfit(ged(x, marks, meth, maxIter=maxIter)), paste("Algorithm ", dQuote(k), " works with markers only."))
		else if( !anyRequired(meth) ) 
			checkTrue(isNMFfit(ged(x, nbasis(model), meth, maxIter=maxIter)), paste("Algorithm ", dQuote(k), " works with no input data."))

	})
	
}

test.autoselect <- function(){
	
	# random data
	x <- rmatrix(100,20)
	model <- rnmf(3, x)
	marks <- gMarkerList(3,5)
	
	.check <- function(.msg, res, selected){
		
		msg <- function(...) paste(.msg, '-', ...)
		checkTrue( isNMFfit(res), msg('plain call works'))
		checkIdentical( algorithm(res), selected, msg('selected algorithm is ', algorithm(res), ': should be ', sQuote(selected)))
	}
	
	# from proportions
	.check('From proportions', ged(x, coef(model)), 'csSAM')	
	.check('From proportions + maxIter=1', ged(x, coef(model), maxIter=1), 'csSAM')
	.check('From proportions + maxIter=1 + markers', ged(x, coef(model), data=marks, maxIter=1), 'cs-qprog')
	if( require(RcppOctave) ){
		.check('From proportions + maxIter>1', ged(x, coef(model), maxIter=10), 'DSection')
	}
	# from signatures
	.check('From signatures', ged(x, basis(model)), 'lsfit')
	.check('From signatures + maxIter=1', ged(x, basis(model), maxIter=1), 'lsfit')
	.check('From signatures + maxIter>1', ged(x, basis(model), maxIter=10), 'lsfit')
	# from markers
	.check('From markers', ged(x, marks), 'DSA')
	.check('From markers + maxIter = 1', ged(x, marks, maxIter=1), 'DSA')
	.check('From markers + maxIter > 1', ged(x, marks, maxIter=10), 'ssKL')
	
}

