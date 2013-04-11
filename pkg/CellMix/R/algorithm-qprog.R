# Partial GED by least-squares using quadratic programming
#
# Optimal quadratic programming [Gong et al. (2011)]
# 
# Creation: 18 Jan 2012
###############################################################################

#' @include registry-algorithms.R
NULL

################################################################################
# CONSTRAINTS
################################################################################

#' Building Constraint Matrices
#' 
#' @description
#' These are utility functions to work with constraint matrices, which define 
#' linear equality and inequality constraints in optimization problems.
#'
#' \code{maxConstraint} creates a constraint matrix that specifies that a 
#' given component must be greater than all the others.
#' It is used, in particular, in the deconvolution method from 
#' \cite{Gong2011} to enforce which is\dQuote{qprog} 
#' 
#' @details 
#' \code{maxConstraint} creates a constraint matrix that specifies that a 
#' component \eqn{k} must be greater than all the others: 
#' 
#' \deqn{x_k \geq \alpha x_i }{x_k >= alpha * x_i},
#' 
#' for all \eqn{1\leq i\leq n}{1<= i <= n}, where \eqn{\alpha}{alpha} is 
#' a fixed parameter is the minimum ratio between the \eqn{k}-th and 
#' all the others components, i.e. it controls the degree of differences 
#' between these components. 
#' 
#' @param n number of components.
#' @param i_max index of the component that must be greater that all others.
#' @param alpha minimum ratio between the \code{k}-th and all the other 
#' components.
#' 
#' @export
#' @rdname constraints
#' @keywords internal
#' 
#' @examples
#' 
#' maxConstraint(5, 3)
#' maxConstraint(5, 3, alpha=3)
#'  
#' 
maxConstraint <- function(n, i_max=1, alpha=1){	
	x <- diag(-alpha, n)
	x[, i_max] <- 1
	x
	x[-i_max, , drop=FALSE]
}

## #' The function \code{toConstraints} converts an object Building Constraint Matrices
## #' 
## #' 
## #' @export
## #' @rdname constraints
## setGeneric('toConstraints', function(object, ...) standardGeneric('toConstraints'))
## #' Creates a constraint matrix from an index vector
## setMethod('toConstraints', 'MarkerList'
##     , function(object, ...){
## 
##     }
## )

################################################################################
# ALGORITHMS
################################################################################

#' Multivariate Least Squares with Equalities and Inequalities
#' 
#' This function provides a multivariate version of function \code{\link{lsei}} 
#' from the \code{\link{limSolve}} package.
#' It solves the following optimisation problem for each column b of B
#' 
#' \deqn{\min_x ||Ax - b||^2}{min ||Ax - b||^2}
#' 
#' subject to the equality constraints: E x = f 
#' 
#' and/or inequality constraints: G x >= h  
#' 
#' @details
#' It applies \code{lsei} on each column of a right-hand side matrix, and 
#' automatically transpose the input and result matrices back and forth if 
#' needed.
#' 
#' @param A numeric matrix containing the coefficients of the quadratic 
#' function to be minimised, \eqn{||Ax-B||^2}; if the columns of \code{A}
#' have a names attribute, they will be used to label the output.
#' @param B numeric matrix containing the right-hand side of the 
#' quadratic function to be minimised.
#' @param fulloutput if \code{TRUE}, also returns the covariance matrix 
#' of the solution and the rank of the equality constraints -- only used 
#' if \code{type = 1} (see \code{\link{lsei}}).
#' @param verbose logical that toggles verbose messages
#' @param ... extra arguments all passed to \code{\link{lsei}}.
#' 
#' @return
#' If \code{fulloutput=TRUE}, a list whose elements are the 
#' results returned by \code{lsei} on each column of \code{B}, 
#' otherwise a list with the following elements:
#' \item{X}{a matrix of dimension \code{ncol(B) x ncol(A)} that 
#' contain the estimated solution of the least-square problem.}
#' \item{IsError}{a logical that indicates if an error occured}
#' \item{type}{the string "lsei", such that how the solution was obtained
#' can be traced} 
#' \item{solutionNorm}{a vector, whose elements are the sum of absolute 
#' values of residuals of equalities and violated inequalities, 
#' as computed for each column of \code{B}.}
#' \item{solutionNorm}{a vector, whose elements are the values of the minimised 
#' quadratic function at the solution, i.e. the value of ||Ax-b||^2, 
#' as computed for each column of \code{B}.}
#' 
#' @import limSolve
#' @import quadprog
#' @export
#' @examples
#' 
#' # random left-hand quadratic coefficients 
#' A <- rmatrix(20, 3)
#' r <- ncol(A) 
#' # random right-hand target matrix
#' B <- rmatrix(20, 10)
#' # constrained least-squares: nonnegative and sum up to one 
#' res <- mlsei(A,B, E=matrix(1, 1, r), F=1, G=diag(1, r), H=rep(0,r))
#' 
mlsei <- function(A, B, ..., fulloutput=FALSE, verbose=TRUE){
	
	library(limSolve)
	residualNorm <- solutionNorm <- numeric()
	IsError <- logical()
	res <- apply(B, 2L, FUN=function(b, ...){
				q <- limSolve::lsei(A=A, B=b, ..., fulloutput=fulloutput, verbose=verbose)
				
				# stack up result fields
				IsError <<- c(IsError, q$IsError)
				residualNorm <<- c(residualNorm, q$residualNorm)
				solutionNorm <<- c(solutionNorm, q$solutionNorm)
				
				if( !fulloutput ) q$X
				else q
			}, ...)
	
	# return plain full output if requested
	if( fulloutput ) return(res)
	
	# return result as a list
	list(X=res, IsError=IsError, type='lsei'
			, residualNorm=residualNorm, solutionNorm=solutionNorm)	
}

#' Proportion or Signature Estimation by Quadratic Programming
#' 
#' This is the workhorse function for the ged algorithms \sQuote{qprog}, 
#' \sQuote{cs-qprog}, and \sQuote{DSA}.
#' 
#' The function implements deconvolution methods that use quadratic programming 
#' techniques to estimate proportions from signatures 
#' (algorithm \code{\link[=qprog-ged]{qprog}} from \cite{Gong2011}), 
#' or both signatures and proportions from marker genes 
#' (algorithm \code{\link[=DSA-ged]{DSA}} from \cite{Zhong2013}).
#' 
#' @param X target matrix
#' @param seed initial NMF model as an \code{\link{NMF}} object.
#' @param data markers as \code{\link{MarkerList}} object.
#' @param exact logical that specifies if one should impose a 
#' sum-to-one (\code{TRUE}) or sum-to-less-than-one (\code{FALSE}) 
#' constraint on the proportions.
#' @param log indicates if the data should be converted to 
#' linear scale (not log-trasnformed).
#' If \code{TRUE} the data is exponentialised (using base 2).
#' If \code{FALSE} the data is left unchanged (the detected log scale is
#' displayed in verbose mode).
#' If a number, then it is used as the base to exponentialise the data.
#' @param ... extra parameters passed to \code{\link{mlsei}}.
#' 
#' @return an \code{\link{NMF}} object.
#' 
#' @keywords internal
.qprog <- function(X, seed, data=NULL, exact=TRUE, log=NULL, ...){
	
	# number of types
	r <- nbasis(seed)
	
	if( hasBasis(seed) && !hasCoef(seed) ){# estimate H from min || Wh - V ||^2
		
		vmessage("Estimating mixture coefficients from basis matrix [method: qprog - "
				, if( exact ) 'exact' else 'inexact',"]")
		## Use constraints for proportions
		# nonnegative: I * P >= 0
		# sum up to one: 1_r * P = 1 (exact) 
		# or sum up to <= one: - 1_r * P >= -1 (inexact)
		res <- 
			if( exact ){
				mlsei(A=basis(seed), B=X
					, E=matrix(1, 1, r), F=1 # equality constraints
					, G=diag(1, r), H=rep(0,r), ...) # inequality constraints
			}
			else{
				mlsei(A=basis(seed), B=X
					, E=NULL, F=NULL # equality constraints
					, G=rbind(diag(1, r), rep(-1, r)), H=c(rep(0,r), -1) # inequality constraints
					, type=2, ...)
			}
		
		coef(seed) <- res$X
		
	}else if( !hasBasis(seed) && hasCoef(seed) ){ # estimate t(W) from min || t(H)w - t(V) ||^2
		vmessage("Estimating basis matrix from mixture coefficients [qprog]")
		
		# will fit: Ax - X'
		lmessage(2L, "Building data range constraint matrices: ", appendLF=FALSE)
		A <- t(coef(seed))
		tX <- t(X)
		# common constraints: x >= min (>= 0) and -x >= -max 
		G.range <- rbind(diag(1,r), diag(-1,r))
		H.range <- c(rep(max(min(X), 0), r), rep(-max(X), r))
		lmessage(2L, round.pretty(H.range[1L]), " <= x <= ", round.pretty(-tail(H.range, 1L)) )
		## Use marker constraints if available
		if( !is.null(data) ){
			if( !isMarkerList(data) )
				stop("qprog - Invalid argument `data`: must be a MarkerList object.")
			if( length(data) != r )
				stop("qprog - Incompatible marker list: length [",length(data),"] is not equal to number of basis component [", r, "]")
			# extract indexes to build constraint matrix
			im <- matchIndex(data, X)
			vmessage("Using markers to build constraint matrices for `mlsei` [", nmark(im), " markers]")
			
			res <- NA_Matrix(r, nrow(X))
			# estimate constrained blocks
			sapply(seq_along(im), function(i){
						# build constraint matrix for the i-th component
						idx <- im[[i]]
						G <- maxConstraint(r, i)
						H <- rep(0, r-1)
						# fit partial target
						q <- mlsei(A, tX[, idx,drop=FALSE]
								, E=NULL, F=NULL # no equality constraints
								, G=rbind(G.range, G), H=c(H.range, H) # inequality constraints: range + marker expression pattern 
								, type=2)
						res[,idx] <<- q$X 
					})
			# estimate "free" blocks
			uidx <- unlist(im)
			q <- mlsei(A, tX[, -uidx,drop=FALSE]
					, E=NULL, F=NULL # no equality constraints
					, G=G.range, H=H.range # inequality constraints: range
					, type=2)
			res[,-uidx] <- q$X 
			
		}else{
			vmessage("Not using any marker constraints")
			
			# estimate
			q <- mlsei(A, tX, E=NULL, F=NULL, G=G.range, H=H.range, type=2)
			res <- q$X
		}
		
		# check for missing values
		if( anyMissing(res) )
			warning("qprog - The estimated basis matrix still contains missing values [", nMissing(res), ']')
		
		# update basis matrix 
		basis(seed) <- t(res)
		
	}else if( isMarkerList(data) ){
		vmessage("Estimating basis and mixture coefficients matrices from marker features [DSA]")
		if( !isMarkerList(data) )
			stop("qprog - Invalid argument `data`: must be a MarkerList object.")
		if( length(data) != r )
			stop("qprog - Incompatible marker list: length [",length(data),"] is not equal to number of basis component [", r, "]")
		# estimate average marker cell-specific expression for each cell type
		# extract indexes to build constraint matrix
		im <- matchIndex(data, X)
		nm <- nmark(im, TRUE)
		vmessage("Using ", sum(nm), "/", nmark(data), " markers to estimate cell proportions: \n"
				, paste0(" ", capture.output(print(nm)), collapse="\n"))
		
		# IMPORTANT: this method works in linear scale only!!!
		vmessage('Checking data scale ... ', appendLF=FALSE)
		if( is_logscale(X) ){
			vmessage('NOTE [log]')
			lbase <- 2
			if( isNumber(log) ){
				lbase <- log
				log <- TRUE
			}
			if( is.null(log) || isTRUE(log) ){				
				vmessage('Converting data to linear scale ... ', appendLF=FALSE)
				X <- expb(X, lbase)
				vmessage('OK [base: ', lbase,']')
			}
		}else{
			vmessage('OK [linear]')
		}
		
		# compute mean expression of each marker set
		p <- .DSAproportions(X, im)
#		if( exact ){ # force exact proportions (rescale sum up to one)
#			p <- scoef(p)
#		}
		# store proportions
		coef(seed) <- p
		##
		# estimate signatures using another call to qprog
		# FIRST: force removing basis matrix to avoid infinite recursive call
		if( hasBasis(seed) )
			basis(seed) <- NA_Matrix(0, nbasis(seed))
		seed <- .qprog(X, seed, data=NULL, ...)
		##
		
	}else
		stop("qprog - Only partial fitting is currently supported (i.e. either the basis or coef matrix is estimated from the other).")
	
	# return updated seed
	seed
}

# Registration of NMF method 'qprog'
nmfAlgorithm.qprog <- setNMFMethod('qprog', .qprog
		, objective='euclidean'
		, mixed=TRUE
		, overwrite=TRUE
)

#' Partial Gene Expression Deconvolution by Quadratic Programming
#' 
#' The GED method \dQuote{qprog} estimates cell/tissue proportions 
#' given a known set of cell/tissue-specific expression signatures, 
#' using quadratic programming techniques as proposed by \cite{Gong2011}.
#' 
#' @inheritParams .qprog
#' @aliases qprog-ged
#' @examples 
#' 
#' ## Example on dummy/random data
#' x <- rmix(3, 100, 20)
#' # true cell signatures
#' s <- basis(x)
#' 
#' # deconvolve using quadratic programming
#' res <- ged(x, s, 'qprog')
#' profplot(x, res)
#' # signatures are not updated
#' identical(basis(res), s)
#' \dontshow{ 
#' 	stopifnot( identical(basis(res), s) )
#' 	stopifnot( nmf.equal(res, ged(x, s, 'qprog')) ) 
#' }
#' 
#'
gedAlgorithm.qprog <- setGEDMethod(key='qprog'
		, description = "Estimates proportions from known expression signatures using quadratic programming"
		, algorithm = 'qprog'
		, reqBasis = TRUE
		, reqCoef= FALSE
		, reqMarker = FALSE
		, maxIter = 1L
		, cite = "Gong2011"
)

#' The GED method \dQuote{cs-qprog} is an \strong{experimental} method that 
#' estimates cell/tissue-specific expression signatures given 
#' cell/tissue proportions, using quadratic programming techniques, 
#' which can optionally impose marker constraints on the signatures.
#'  
#' @aliases cs-qprog-ged
#' @rdname gedAlgorithm.qprog
#' @examples 
#' 
#' ### CS-QPROG
#' ## Example on dummy/random data
#' # random target matrix
#' x <- rmix(3, 100, 20)
#' # true cell proportions
#' p <- coef(x)
#' 
#' # deconvolve using quadratic programming
#' res <- ged(x, p, 'cs-qprog')
#' profplot(basis(x), basis(res))
#' 
#' # signatures are not updated
#' identical(coef(res), p)
#' \dontshow{ 
#' 	stopifnot( identical(coef(res), p) )
#' 	stopifnot( nmf.equal(res, ged(x, p, 'cs-qprog')) ) 
#' }
#' 
#'
gedAlgorithm.cs_qprog <- setGEDMethod(key='cs-qprog'
		, description = "Estimates constrained cell-specific signatures from proportions using quadratic programming [experimental]"
		, algorithm = 'qprog'
		, reqBasis = FALSE
		, reqCoef= TRUE
		, reqMarker = FALSE
		, maxIter = 1L
)
