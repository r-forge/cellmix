# Unit tests for class MarkerList
# 
# Author: Renaud Gaujoux
# Creation: 14 Jan 2012
###############################################################################

src.data <- function(rfun=gmarkers){
	
	l_case <- list()
	.case <- function(msg, val){
		l <- list(msg=msg, val=val)
		l_case <<- c(l_case, list(l))
	}
	.case("Marker vector with NO names NO values", rfun(3))
	.case("Marker vector with only names", rfun(3, names=TRUE))
	.case("Marker vector with values only", rfun(3, values=TRUE))
	.case("Marker vector with values and names ", rfun(3, names=TRUE, values=TRUE))
	.case("Marker vector with names and index values", rfun(3, names=TRUE, values=1:18))
	
	l_case
}

#' Test random generation of markers
test.rmarkers <- function(){
	
	.MSG <- NULL
	msg <- function(...) paste(.MSG, '-', ...)
	
	.check <- function(m, n, multfactor=-3){
		checkTrue(!anyMissing(m), msg("result has no missing values"))
		checkEquals( length(m), n*(n+1)/2*abs(multfactor), msg("result has the correct length"))
		
		if( is.character(m) ){
			checkEquals( nlevels(factor(m)), n, msg("result has the correct number of types"))
			checkEquals( summary(factor(m)), 1:n * abs(multfactor), msg("result has the correct breakdown of types")
						, checkNames=FALSE)
		}
		if( !is.null(names(m)) ){
			checkTrue(!anyDuplicated(names(m)), msg("result has no duplicated names"))
			checkTrue(!anyMissing(names(m)), msg("result has no missing names"))
			checkTrue(!any(names(m)==''), msg("result has no empty names"))
		}
	}
	
	# single argument
	.MSG <- "Single argument"
	m <- rmarkers(3)
	checkTrue(is.character(m), msg("result is a character vector"))
	.check(m, 3)
	
	# names = TRUE
	.MSG <- "With names=TRUE"
	m <- rmarkers(3, names=TRUE)
	checkTrue(is.character(m), msg("result is a character vector"))
	checkTrue(!is.null(names(m)), msg("result has names"))
	checkIdentical(length(unique(names(m))), 3L, msg("result correct number of types"))
#	.check(m, 3)
	
	# names a character vector
	.MSG <- "With names a character vector"
	checkException(rmarkers(3, names=letters[1:3]), msg("Error if names not long enough"))
	checkException(rmarkers(3, names=rep(letters[1:3],100)), msg("Error if names not long enough, after removing duplicates"))
	.MSG <- "With names a character vector with exact number of names"
	m <- rmarkers(3, names=letters[1:20])
	checkTrue(is.character(m), msg("result is a character vector"))
	checkTrue(!is.null(names(m)), msg("result has names"))
	checkIdentical(length(unique(names(m))), 3L, msg("result correct number of types"))
#	.check(m, 3)
	
	# with values
	.MSG <- "With values a numeric vector"
	v <- runif(3)
	m <- rmarkers(3, values=v)
	checkIdentical(sort(unique(unlist2(m))), sort(v), 'Values are sampled with replacement if not long enough') 
	.MSG <- "With values TRUE"
	m <- rmarkers(3, values=TRUE)
	checkTrue(is.list(m), msg("result is a list"))
	.check(unlist(m), 3)
	
	
}

test.stack <- function(){
	
	.MSG <- NULL
	msg <- function(...) paste(.MSG, '-', ...)
	
	.checkRes <- function(x, ref, ...){
		checkTrue(is.data.frame(x), msg('result is a data.frame'))
		checkEquals(colnames(x)
			, c('values', 'ind', 'names')
			, msg("result has 3 columns with names 'values', 'ind' and 'names'") )
		
		checkEquals(x$names, setNames(marknames(ref), NULL), msg("Column `names` contains the marker names"))
		if( hasValues(ref) )
			checkEquals(x$values, unlist(ref, use.names=FALSE), msg("Column `values` contains the values"))
		else 
			checkTrue(all(is.na(x$values)), msg("Column `values` contains only NAs"))
	}
	
	l_case <- list()
	.case <- function(msg, val){
		l <- list(msg=msg, val=val)
		l_case <<- c(l_case, list(l))
	}
	.case("List with NO names NO values", rmarkers(3))
	.case("List with only names", rmarkers(3, names=TRUE))
	.case("List with values only", rmarkers(3, values=TRUE))
	.case("List with values and names ", rmarkers(3, names=TRUE, values=TRUE))
	lapply(l_case, function(x){
		.MSG <<- x$msg
		ml <- MarkerList(x$val)
		.checkRes(stack(ml), ml)
	})
	
}

#test.flatten <- function(){
#	
#	.MSG <- NULL
#	msg <- function(...) paste(.MSG, '-', ...)
#	input <- src.data()
#	lapply(input, function(x){
#		.MSG <<- x$msg
#		checkIdentical(x$val, flatten(MarkerList(x$val))
#				, msg("Compounding flatten and MarkerList leads to identity on [", str_out(x$val, use.names=TRUE), "]"))
#		
#	})
#	
#}

#' Tests for constructor function MarkerList
test.MarkerList <- function(){
	
	.MSG <- NULL
	msg <- function(...) paste(.MSG, '-', ...)
	
	# input
	input <- src.data()
	lapply(input, function(x){
		.MSG <<- x$msg
		checkTrue(is(ml <- MarkerList(x$val), 'MarkerList')
		, msg("Constructor MarkerList on [", str_out(x$val, use.names=TRUE), "] returns a MarkerList object."))
		
	})
	
} 

#' Tests method `[`
test.brackets <- function(){
	
	set.seed(1)
	# define reference data	
	x <- lapply(seq(5,30,by=5), function(n) setNames(runif(n), paste('name', n, 1:n,sep='_')))
	names(x) <- toupper(letters[1:length(x)])
	m <- MarkerList(x)
	
	.local <- function(msg, ...){
		
		.msg <- function(...) paste(msg, ': ', ..., sep='')
		
		ms <- m[...]
		x <- x[...]
		checkTrue( is(ms, 'MarkerList'), .msg('Result is a MarkerList object'))
		checkIdentical(ms@.Data, x, msg("Data is correctly subset"))
		checkIdentical(names(ms), names(x), msg("Names are correctly subset"))
		sapply(slotNames(m), 
			function(s) checkIdentical(slot(ms, s), slot(m, s), .msg("Slot '", s, "' is unchanged"))
		)
	}
	
	i_case <- data.frame(msg=NA, val=NA, stringsAsFactors=FALSE)
	.case <- function(msg, val) i_case <<- rbind(i_case, c(msg=msg, val=val))
	# single integer
	.case('i=single integer [1]', 1)
	.case('i=single integer [3]', 3)
	.case('i=single integer [end]', length(x))
	# integer vector
	.case('i=ordered integer vector from 1', 1:3)
	.case('i=ordered integer vector from >1', 3:5)
	.case('i=unordered integer vector', c(2,4,3,1))
	.case('i=unordered integer vector with duplicates', c(2,4,2,3))
	# TODO: character strings
  # TODO: lists
}

checkIsMarkerList <- function(x, ...) checkTrue(isMarkerList(x), ...)

test.subset <- function(){
	
	m <- MarkerList(letters[1:20], names=gl(4, 5, labels=paste0('Type', 1:4)))
	sub <- letters[20:15]
	ref_ul <- c(letters[15], letters[16:20])
	oref_ul <- c(letters[15], letters[20:16]) 
			
	checkIsMarkerList(m2 <- subset(m, sub), 'Subset by character vector returns a MarkerList object')
	checkIdentical(names(m2), names(m), 'Subset does not change element names/order')
	checkIdentical(as.numeric(summary(m2)), c(0, 0, 1, 5), 'Subset gives correct breakdown')
	checkEquals(unlist2(m2), ref_ul, checkNames=FALSE, 'Subset is correct')
	
}