# Utilities for NMF constraints
# 
# Author: Renaud Gaujoux
# Creation: 15 Nov 2011
###############################################################################


#' Checking Linear Constraints
#' 
#' Checks if a given matrix or vector satisfies a set of linear constraints.
#' The constraints are specified by a left-hand matrix \eqn{A} and a right-hand 
#' vector \eqn{b} such that a matrix or vector \eqn{x} satisfies the constraints 
#' if an only if \deqn{A \cdot x \geq b}{A x >= b}.
#' 
#' @param x matrix or vector to check.
#' @param A left-hand matrix
#' @param b right-hand vector
#' @param margin which dimension should be checked (1L: rows or 2L: columns).
#' @param op comparison operator to use (\code{'>=', '<=', '='}), indicating the 
#' type of constraints: inequality (greater than or lower than) or equality.
#' 
#' @return a logical vector of length \code{nrow(A)}.
#' 
#' @export 
#' 
checkConstraints <- function(x, A, b, margin=2L, op = c('>=', '<=', '==')){
	
	if( length(A) == 0 ) return(TRUE)
	
	op <- match.arg(op)
	op <- 
	if( op == '==' ) function(a, b) all.equal(a, b, check=FALSE)		
	else getFunction(op)
	val <- if( margin == 2L || is.vector(x) ) A %*% x else tcrossprod(A, x)
	apply(val, 2L, function(v) op(v, b))
	
}
