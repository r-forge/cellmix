# GED algorithm common development utilities 
# 
# Author: Renaud Gaujoux
# Created: Mar 26, 2013
###############################################################################

#' @include registry-algorithms.R
NULL

#' @importFrom corpcor pseudoinverse 
.robust_pseudoinverse <- function(x, ...){
	res <- pseudoinverse(x, ...)
	res[is.na(res)] <- 0.001
	res
} 

.lsqnonneg <- function(x, y, verbose=FALSE, check=TRUE){
	
	# call the internal function
	if( check && rcond(crossprod(x)) < .Machine$double.eps ){
#		message('x')
#		print(head(x))
#		message('pseudo')
#		print(head(.robust_pseudoinverse(xtx)))
		eps <- .Machine$double.eps
		w <- which(x < eps)
		if( length(w) ) x[w] <- 10^-6 # runif(length(w), 0, 10^-6)
	}
	res <- NMF::.fcnnls(x, y, verbose = verbose, pseudo = .robust_pseudoinverse
	#, eps = .Machine$double.eps * 10 * norm(x, '1') * length(x) 
	)
	
	# put dimnames back in if necessary
	if( is.null(rownames(res$coef)) ) rownames(res$coef) <- colnames(x)
	if( is.null(colnames(res$coef)) ) colnames(res$coef) <- colnames(y)
	
	# wrap up the result
	list(x=res$coef, resnorm = norm(x %*% res$coef - y, '2'))
	
}

# complete fast seeding of cell-specific signatures using marker indexes
# - qprog for enforcing marker patterns
# - fcnnls for the unconstrained features
csfit_markers <- function(x, p, markers, ratio=2, G=NULL, H=NULL, verbose=FALSE){
	
	if( verbose ){
		message("Computing constrained cell-specific signatures using qprog/fcnnls method"
				, " (cs-fold: ", ratio, ") ... ", appendLF=FALSE)
	}
	
	r <- nrow(p)
	stopifnot( length(markers) == r )
	A <- t(p)
	X <- x
	tX <- t(X)
	im <- markers
	
	# common constraints: x >= min (>= 0) and -x >= -max
	G.range <- G
	H.range <- H
	if( is.null(G) ){
		G.range <- rbind(diag(1,r), diag(-1,r))
		H.range <- c(rep(0, r), rep(-max(X), r))
	}
	
	# init result
	res <- NA_Matrix(r, nrow(X))
	# put dimnames back in if necessary
	colnames(res) <- rownames(x)
	rownames(res) <- names(markers)
	
	# estimate constrained blocks
	sapply(seq_along(im), function(i){
				# build constraint matrix for the i-th component
				idx <- im[[i]]
				G <- maxConstraint(r, i, alpha=ratio)
				H <- rep(0, r-1)
				# fit partial target
				q <- mlsei(A, tX[, idx,drop=FALSE]
						, E=NULL, F=NULL # no equality constraints
						, G=rbind(G.range, G), H=c(H.range, H) # inequality constraints: range + marker expression pattern 
						, type=2)
				res[,idx] <<- q$X 
			})
	
	res <- t(res)
	
	# estimate "free" blocks - if any -- using fcnnls
	if( length(uidx <- unique(unlist(im))) < nrow(x) ){				
		q <- .lsqnonneg(A, tX[, -uidx, drop=FALSE], check=FALSE)
		res[-uidx, ] <- q$x
#		res <- .lsqnonneg(A, tX[,,drop=FALSE], check=FALSE)$x
	}
	
	if( verbose ){
		message("OK [", length(uidx), " constrained - ", nrow(x) - length(uidx), " free]")
	}
	
	# return estimated cs-signatures
	res
}
