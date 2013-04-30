# Unit tests for matrix stats
# 
# Author: Renaud Gaujoux
# Created: Apr 30, 2013
###############################################################################



test.applyBy <- function(){
	
	# random data matrix
	x <- rmatrix(12, 6)
	# groups of columns
	fc <- gl(2, 3)
	fr <- gl(3, 4)
	
	
	.check <- function(AFUN, x, margin, gr, resclass = class(x), dims, resObject=NULL, msg=NULL){
		
		.msg <- function(...){
			paste0(if( !is.null(msg) ) paste0(msg, ' - '), "Function '", AFUN,"' + input '", class(x), "' + margin=", margin, ': ', ...)
		}
		FUN <- match.fun(AFUN)
		
		check_alt <- FALSE
		b <- if( is.null(resObject) ){
			check_alt <- TRUE
			applyBy(x, gr, margin, FUN) 
		}else resObject
		checkTrue(is(b, resclass), .msg("result is a '", resclass, "'"))
		checkIdentical( setNames(dim(b)[c(margin, 3L-margin)], NULL), setNames(c(dim(x)[margin], nlevels(gr)), NULL), .msg('output dimensions are correct'))
		checkIdentical(dimnames(b)[[margin]], dimnames(x)[[margin]], .msg("margin names are correct (unchanged)"))
		checkIdentical(dimnames(b)[[3L - margin]], levels(gr), .msg("by margin names are correct (levels of by factor)"))
		
		# check alternative call
		if( check_alt ){
			altcall <- if( margin == 1L ) 'rowApplyBy' else 'colApplyBy' 
			balt <- do.call(altcall, list(x, gr, FUN))
			checkIdentical(exprs(b), exprs(balt), .msg('result is identical to alternative no-margin call "', altcall, '"'))
		}
		
		if( is(x, 'ExpressionSet') ){
			checkIdentical(annotation(b), annotation(x), .msg('Annotation are conserved'))
			if( margin == 1L ){				
				checkIdentical(fData(b), fData(x), .msg('Feature data are conserved'))
				checkIdentical(nrow(pData(b)), nlevels(gr), .msg('Sample data are collapsed'))
			}else{
				checkIdentical(nrow(fData(b)), nlevels(gr), .msg('Feature data are collapsed'))
				checkIdentical(pData(b), pData(x), .msg('Sample data are conserved'))
			}
		}
		
		# return result
		b
	}
	
	.checkAll <- function(x, ...){
		.check('rowSums', x, 1L, fc, ...)
		.check('colSums', x, 2L, fr, ...)
	}
	
	# matrix
	.checkAll(x)
	# expression set	
	x <- ExpressionSet(x, annotation='abcd.db')
	.check('rowSums', x, 1L, fc)
	.check('colSums', x, 2L, fr)

	## annotations are conserved/collapsed
	pData(x) <- data.frame(Group=rev(fc), Sample=letters[1:ncol(x)])
	fData(x) <- data.frame(ENTREZID=rev(fr), Gene=letters[nrow(x):1])
	y <- .check('rowSums', x, 1L, fc, resObject = applyBy(x, 'Group', 1L, rowSums), msg = 'Using phenotypic variable name')
	checkTrue(all(rownames(pData(y)) == as.character(y$Group)), 'Effect of rowSums on pheno data: Collapsed by variable is consistent with the aggregated sample names')
	y <- .check('colSums', x, 2L, fr, resObject = applyBy(x, 'ENTREZID', 2L, colSums), msg = 'Using feature variable name')
	checkTrue(all(rownames(fData(y)) == fData(y)$ENTREZID), 'Effect of colSums on feature data: Collapsed by variable is consistent with the aggregated feature names')
	TRUE
}