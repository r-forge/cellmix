
#' @include AllGenerics.R 
NULL

#' Class for Gene Expression Deconvolution Benchmark Datasets
#' 
#' @slot pure logical vector that indicates which samples are pure. 
#' 
#' @export
setClass('ExpressionMix', contains=c('ExpressionSet', 'NMFstd')
	, representation = representation(
		pure = 'logical'
	)
	, validity = function(object){
		
		# the mixture dimensions must be compatible with the expression data
		dm <- selectMethod('dim', 'NMF')(object)[1:2]
		de <- selectMethod('dim', 'ExpressionSet')(object)		
		if( hasBasis(object) && de[1] != dm[1] )
			return(paste("Incompatible number of rows: mixture [", dm[1],"]"
				, " vs. expression [", de[1],"]", sep=''))
		if( hasCoef(object) && de[2] != dm[2] )
			return(paste("Incompatible number of columns: mixture [", dm[2],"]"
					, " vs. expression [", de[2],"]", sep=''))

		# the mixture proportions must sum up to 1
		if( nbasis(object) > 0 && ncol(object) > 0 ){
			co <- coef(object)
			p <- colSums(co)
			# check for missing values => throw warning if any
			if( length(nna <- which(is.na(p))) > 0 )
				warning("Missing proportions for ", length(nna), " sample(s): "
								, paste(head(sampleNames(object)[nna], 3), collapse=', '), ", ... [check coef(object)]")
			if( length(nok <- which(abs(p-1) > 0.01)) > 0 )
				return(paste("Invalid proportions for ", length(nok), " sample(s) [should sum up to one +- 1%]: "
					, paste(head(sampleNames(object)[nok], 3), head(p[nok], 3), sep='=', collapse=", "), ", ...", sep=''))	
		}
		
		# the logical vector `pure` must have the correct length
		if( length(object@pure)>0 && length(object@pure) != ncol(object) )
			return(paste("Invalid slot `pure`: length [", length(object@pure)
					,"] must match the number of samples [", ncol(object),"]", sep=''))
			
		TRUE
	}
)

#' Initialize method for ExpressionMix object
setMethod('initialize', 'ExpressionMix', 
	function(.Object, ...){
				
		# initialize the ExpressonSet object
		.Object <- callNextMethod(.Object, ...)
		
		# add version information
		classVersion(.Object)['ExpressionMix'] <- '1.0.0'
		
		# synchronize dimnames: use the ExpressionSet dimnames as they are always there
		dimnames(.Object) <- list(featureNames(.Object)
								, sampleNames(.Object)
								, basisnames(.Object))
		
		.Object
})


#' Factory Method for ExpressionMix Objects
#' 
#' This is a generic function that creates instances of 
#' \code{\linkS4class{ExpressionMix}} objects from common input data.
#' 
#' @param object main input data. See possible values in section \emph{Methods}.
#' @param ... extra arguments to allow extension.
#' See each method's description for more details.
#' 
#' @return an \code{\linkS4class{ExpressionMix}} object
#' 
#' @export
#' @inline
setGeneric('ExpressionMix', function(object, ...) standardGeneric('ExpressionMix'))


#' Show Data Available in an Object
#' 
#' Shows which data is available in an object.
#' For \code{\linkS4class{ExpressionMix}} objects, the these are amongst 
#' expression, cell signature and/or cell mixing proportion data.
#' 
#' @param object object for which on wants to see the available data
#' @param ... extra arguments for future extension  
#'
#' @inline
#' @export  
setGeneric('showData', function(object, ...) standardGeneric('showData'))


#' Extracting Expression Data
#' 
#' Extracts expression data (usually as \code{ExpressionSet} objects) from 
#' other objects.
#' 
#' @param object object from which to extract the expression data
#' @param ... extra argument to allow extension.
#' See each method's description for more details,.
#' 
#' @inline
#' @export
setGeneric('eset', function(object, ...) standardGeneric('eset') )


#' Extracting Mixture Data
#' 
#' The S4 generic function \code{mixData} extract mixture data 
#' from an object, i.e. information on its underlying constituents/components.
#' For example, for gene expression deconvolution, the \pkg{CellMix} stores such data
#' as an \code{\linkS4class{NMFstd}} model object.
#' 
#' @param object object from which to extract the mixture data
#' @param ... extra argument to allow extension.
#' See each method's description for more details,.
#' 
#' @inline
#' @export
setGeneric('mixData', function(object, ...) standardGeneric('mixData') )
#' The function \code{mixData<-} assign mixture data onto objects.
#' 
#' @param value replacement value.
#' In the case of \code{\link{ExpressionMix}} objects, this will be an 
#' object of class \code{\linkS4class{NMFstd}}.  
#'
#' @rdname mixData 
#' @inline
#' @export
setGeneric('mixData<-', function(object, value) standardGeneric('mixData<-') )


## AUTO-COMPLETION 
setGeneric('.DollarNames', package='utils')

eset_completion_matches <- function(object, pattern = ""){
	grep(pattern, varLabels(object), value=TRUE)	
}

#' @S3method .DollarNames ExpressionMix
.DollarNames.ExpressionMix <- function(x, pattern = "") eset_completion_matches(x, pattern)

#' Auto-completion for \code{\linkS4class{ExpressionMix}} objects
#' @export
setMethod('.DollarNames', 'ExpressionMix', .DollarNames.ExpressionMix)
