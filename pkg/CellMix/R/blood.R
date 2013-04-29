# Blood-related functions and data
#
# Author: Renaud Gaujoux
# Creation: 15 Jan 2013

#' @include ged.R
NULL

#' Blood Sample Deconvolution
#' 
#' @description
#' The functions described here are dedicated to gene expression deconvolution
#' of blood samples (i.e. whole blood or PBMCs).
#' 
#' \code{gedBlood} uses the methodology defined by \cite{Abbas2009}, which uses a 
#' fixed set of 17 cell type-specific signatures to estimate cell proportions in 
#' blood samples.
#' Each signature corresponds to a white blood cell in resting or activated state
#' (See section Details).
#' 
#' @details
#' 
#' The signatures used by \code{gedBlood} were designed by \cite{Abbas2009} to optimise 
#' their deconvolution power.
#' They are available in the \pkg{CellMix} as dataset \code{\link{Abbas}}.
#' 
#' \code{gedBlood} is essentially a shortcut for \code{gedProportions(object, Abbas, ...)}, 
#' see \code{\link{gedProportions}} for details on the possible arguments.
#' 
#' @param object target data, specified in any format supported by \code{\link{ged}}.
#' For \code{asCBC}, an object for with suitable \code{asCBC} method defined.
#' @param ... extra arguments passed to \code{\link{gedProportions}}.
#' @inheritParams gedProportions
#' 
#' @seealso \code{\link{gedProportions}}, \code{\link{ged}}
#' 
#' @rdname blood
#' @export
#' @examples
#' 
#' if( !isCHECK() ){
#' 
#' # load kidney transplant data (Shen-Orr et al. (2010))
#' eset <- ExpressionMix('GSE20300')
#' 
#' # true CBC proportions are stored as -- mixture -- coefficients
#' cbc <- coef(eset)
#' head(cbc)
#' profplot(cbc, lab='')
#' 
#' # devonvolve using basis signature matrix from Abbas et al. (2009)
#' res <- gedBlood(eset, verbose=TRUE)
#' # estimated proportions are stored as -- mixture -- coefficients
#' p <- coef(res) 
#' str(p)
#' # the Abbas basis matrix includes detailed immune cell types
#' rownames(p)
#' # or: basisnames(res)
#' 
#' # aggregate into CBC
#' ecbc <- asCBC(res)
#' # plot Estimated vs. Measured CBC 
#' profplot(cbc, ecbc)
#' 
#' }
#' 
gedBlood <- function(object, ..., normalize=TRUE, verbose=FALSE){
	
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	verbose <- lverbose()
	
	# load Abbas signatures
	vmessage("Loading basis signature from Abbas et al. (2009) ... ", appendLF=FALSE)
	abbas <- packageData('Abbas')
	vmessage('OK [', nrow(abbas), ' features x ', ncol(abbas), ' cell types]')
	if( verbose ) logger.options(autoindent=FALSE) 
	# fit proportions
	gedProportions(object, abbas, ..., normalize = normalize, verbose = verbose)
}

# mapper to translate cell type names, e.g., to aggregate blood cell subsets into CBC subsets
.charmap <- function(x, map){
	map <- setNames(as.character(map), names(map))
	i <- match(tolower(x), tolower(names(map)))
	ok <- !is.na(i)
	i[ok] <- as.character(map[i[ok]])
	as.character(i)
}

cbcMapper <- list(
	Abbas = function(x){
		Abbas <- ldata('Abbas')
		.charmap(x, setNames(Abbas$WBType, sampleNames(Abbas)))
	}
	, DMap = list(Lymphocytes=c('NKT', 'CD4+ Central Memory', 'CD4+ Effector Memory'
								, 'CD8+ Central Memory', 'CD8+ Effector Memory', 'CD8+ Effector Memory RA'))
	, HaemAtlas = c(`B-CD19`='Lymphocytes'
			, Erythroblast='Precursor', `Granulocyte-CD66b`='Neutrophils'
			, Megakaryocyte='Precursor', `Monocyte-CD14`='Monocytes', `NK-CD56`='Lymphocytes'
			, `T-CD4`='Lymphocytes', `T-CD8`='Lymphocytes')
	, regexpr = function(x){
		# define patterns
		patterns <- setNames(sub("s$", "", names(refCBC)), names(refCBC))
		patterns <- c(patterns, Lymphocytes='((^)|( ))((B)|(T[ch]?)|(NK)|(natural[- ._]killer)|(NKT))((( )|($))|([- ._]?cell))'
					, Granulocytes='granulocyte')
		patterns <- c(patterns, c(Lymphocytes='^CD[48][^0-9]'))
			
		# apply patterns sequentially
		res <- as.character(rep(NA, length(x)))
		for(ip in seq_along(patterns)){
			if( !length(i <- which(is.na(res))) ) break;
			w <- grepl(patterns[ip], x[i], ignore.case=TRUE)
			res[i[w]] <- names(patterns[ip]) 
		}
		# return result
		res
	}
)

#' \code{asCBC} is a S4 generic function that converts detailed -- high resolution -- 
#' proportions of blood cells into Complete Blood Count (CBC) proportions, 
#' i.e. Monocytes, Basophils, Lymphocytes, Neutrophils and Eosinophils.
#' 
#' \code{asCBC} has methods defined for \code{NMF} models and \code{Markerlist} objects.
#' See each method's description for more details.
#' 
#' Currently \code{asCBC} methods will correctly work only on objects that have cell types 
#' that match exactly names of signatures in the \code{\link{Abbas}} dataset.
#' 
#' @rdname blood
#' @inline
#' @export
setGeneric('asCBC', function(object, ...) standardGeneric('asCBC') )
#' This is the workhorse method that maps immune/blood cell type names 
#' to the CBC cell types: Monocytes, Basophils, Lymphocytes, Neutrophils and Eosinophils.
#' 
#' It returns a factor, whose names are elements of \code{object} and the values are 
#' their corresponding CBC cell type.
#' If \code{drop=FALSE} the result is of the same length as \code{object}, otherwise
#' it only contains elements that could be mappped to a cell type.  
#' 
#' @param drop logical that indicates if elements in \code{object} that cannot be mapped to 
#' a cell type should be removed from the returned mapping.
#' @param quiet logical that indicates that the mapping should be performed quietly.
#' If \code{FALSE}, then an error is thrown if none of the elements can be mapped, 
#' or, if in addition \code{drop=FALSE}, a warning is thrown if only some of the elements
#' could be mapped.  
#'    
setMethod('asCBC', 'character'
	, function(object, drop=FALSE, quiet=FALSE){
		
		# sequentially apply each CBC mapper
		res <- setNames(as.character(rep(NA, length(object))), object)
		for( map in cbcMapper ){
			# stop as soon as all type is mapped
			if( !length(i <- which(is.na(res))) ) break;
			
			# match unmapped type
			ct <- names(res)[i]
			if( is.function(map) ){
				m <- map(ct)
			}else if( is.character(map) ){
				m <- .charmap(ct, map)
			}else if( is.list(map) ){
				map <- unlist2(map)
				m <- .charmap(ct, setNames(names(map), map))
			}else stop("Invalid cell type map [", class(map), ']')
			# update result map
			if( !is.null(m) ) res[i] <- m
		}
		
		cbcTypes <- names(refCBC)
		# check that at least some cell types were matched
		if( all(is.na(res)) ){
			if( !quiet ) stop("asCBC - Could not match any cell type to a CBC type.")
			else if( drop ) return(factor(,levels=cbcTypes))
			else return(factor(res,levels=cbcTypes))
		}
		# check that all cell types were matched
		if( any(bad <- is.na(res)) ){
			# drop unmatched cell types
			if( !quiet ) warning("asCBC - Could not match some cell type(s): ", str_out(unique(names(res[bad])), Inf))
		}
		
		# limit to CBC cell types
		if( drop ) cbcTypes <- cbcTypes[cbcTypes %in% res]
		# return as a factor
		factor(res, levels=cbcTypes)
	}
)
#' The result of gene expression deconvolution performed by \code{\link{ged}} are 
#' stored in \code{\linkS4class{NMFstd}} model objects, 
#' which contain the cell type-specific signatures and/or cell relative 
#' proportions.
#' 
#' This method aggregates, i.e. sums up, the cell proportions and averages 
#' the signatures of cell types from each of the CBC groups that are 
#' available in the data.
#' 
setMethod('asCBC', 'NMF'
	, function(object, drop=TRUE, ...){
		if( is.null(cl <- basisnames(object)) ){
			stop("Cannot compute/aggregate to CBC data: object as no cell type names [basisnames].")
		}
		cbc <- asCBC(cl, ..., drop=drop)
		# sum proportions
		coef(object) <- colSumsBy(coef(object), cbc)
		# average signatures
		basis(object) <- rowMeansBy(basis(object), cbc)
		# TODO: add dropped basis
		# return aggregate object
		object
	}
)
#' This method combines markers of cell types that belong to the same CBC group.  
setMethod('asCBC', 'MarkerList' 
		, function(object, ...){
			cbc <- asCBC(names(object), ...)
			combine(object, cbc)
		}
)

# TODO: add reference to Sarwal paper here

#' \code{refCBC} is a numeric vector that contains average Complete Blood Count proportions (CBC) in healthy persons, 
#' based on empirical studies in healthy patients.
#' It contains proportions for Basophils, Lymphocytes, Eosinophils, Neutrophils and Monocytes.
#' 
#' @rdname blood
#' @export
refCBC <- c(Basophils = 0.5, Lymphocytes = 29.5, Eosinophils = 3, Neutrophils = 57, Monocytes = 10)
refCBC <- refCBC / sum(refCBC)
# Relative proportion based on Sarwal data for CBC. (DATA ORGANIZED IN SAME ORDER AS IN REGEV DATA)

#' \code{gCBC} generates a matrix of average Complete Blood Count proportions (CBC) for a given number of samples.
#' The default proportions are based on empirical studies in healthy patients (see \code{\link{refCBC}}), 
#' and each sample get assigned the same proportions.
#' 
#' @param n number of samples in the generated CBC matrix
#' @param sampleNames names of the samples, recycled or truncated if necessary, to match \code{n}.
#' @param counts CBC data to use instead of the defaults.
#' It must be a numeric vector.
#' 
#' @examples
#' 
#' # default proportions
#' gCBC()
#' 
#' # for 5 samples
#' gCBC(5)
#' 
#' # setting sample names
#' gCBC(5, sampleNames = letters[1:10]) # names are truncated if necessary
#' 
#' @rdname blood
#' @export
gCBC <- function(n=1, sampleNames=NULL, counts=NULL){
	
	# Relative proportion based on Sarwal data for CBC. (DATA ORGANIZED IN SAME ORDER AS IN REGEV DATA)
	if( is.null(counts) )	p <- refCBC
	else counts <- p
		
	# generate matrix
	p <- replicate(n, p)
	# add colnames
	if( !is.null(sampleNames) )
		colnames(p) <- rep(sampleNames, length.out=ncol(p))
	# return scaled proportions
	scoef(p)
}
