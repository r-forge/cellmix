# Unit tests for idtype
# 
# Author: Renaud Gaujoux
# Creation: 24 Apr 2012
###############################################################################


test.idtype <- function(){
	
	.test_id <- function(id, type, rev=FALSE, limit=NULL){
		if( rev ) id <- revmap(id)
		checkTrue(all(idtype(keys(id), limit=limit)==type), str_c("ID type `", type, "` is correctly detected"))
	}
	
	library(illuminaHumanv2.db)
	
	.test_id(illuminaHumanv2ENSEMBL2PROBE, "ENSEMBL")
	.test_id(illuminaHumanv2GO2PROBE, "GOID")
	.test_id(illuminaHumanv2UNIGENE, "UNIGENE", TRUE)
	.test_id(illuminaHumanv2REFSEQ, "REFSEQ", TRUE)
	.test_id(illuminaHumanv2NUID, ".nuID", TRUE)
	
	# gene bank
	checkIdentical(idtype(c('AB123456', 'A12345')), "GENEBANK", "Gene bank nucleotide is ok")
	checkIdentical(idtype(c('ABC12345', 'DEF12345')), "GENEBANK", "Gene bank protein is ok")
	checkIdentical(idtype(c('ABCD12345678', 'EFGH123456789', 'IJKL1234567890')), "GENEBANK", "Gene bank WGS is ok")
	checkIdentical(idtype(c('ABCDE1234567', 'FGHIJ1234567')), "GENEBANK", "Gene bank MGA is ok")
	
#	ids <- list(	
#	.AFFY = c("AFFX-2132asad", "241838_at", "241838_x_at")
#	, ENSEMBL = c("ENS013135", "ENS213216")
#	, ENSEMBLTRANS = c("ENST013135", "ENST213216")
#	, ENSEMBLPROT = c("ENSP013135", "ENSP213216")
#	, ENTREZID = c(13135, 213216)
#	, UNIGENE = c("Hs.123456", "Rn.2346545")
#	, IMAGE = c("IMAGE:131265", "IMAGE:9876")
#	, GOID = c("GO:123", "GO:1")
#	, PFAM = c("PF213", "PF5")
#	)
#	
#	mapply(.test_id, ids, names(ids))
}

test.asGeneIdentifierType <- function(){
	
	fun <- asGeneIdentifierType
	# missing
	checkTrue(is.character(known <- fun()), "missing x: returns character vector of known types")
	# known types, including partial match
	sapply(known, function(x){
		t <- fun(x)
		.msg <- function(...) paste('result for ', x, ': ', ..., sep='')
		checkTrue(is(t, 'GeneIdentifierType'), .msg("is a GeneIdentifier object"))
		checkTrue(!is(t, 'NullIdentifier'), .msg("is NOT a NullIdentifier object"))
		checkIdentical(t, fun(substr(x, 0, 3)), .msg("partial match works"))
		annotation(t) <- 'aaa'
		checkIdentical(fun(t), t, .msg("call on GeneIdentifier object returns identical (and conserve annotation)"))
		checkIdentical(fun(t, annotation='toto'), new(class(t), annotation='toto')
			, .msg("call on GeneIdentifier object + annotation returns same class with annotation changed"))
	})
	#
	checkIdentical(fun('ann.db'), AnnotationIdentifier('ann.db'), 'works with annotation package name')
	checkIdentical(fun('.Affymetrix'), AnnotationIdentifier(), 'works with probe type name (e.g., .Affymetrix)')
	checkIdentical(fun('.Affymet'), AnnotationIdentifier(), 'works with partial match on probe type name (e.g., .Affymet)')
	checkIdentical(fun('.Affymetrix', annotation='toto'), AnnotationIdentifier('toto')
		, 'take annotation into account')
	# extract from data
	g <- GeneSet('a', geneIdType=EntrezIdentifier())
	checkIdentical(fun(g), geneIdType(g), 'extract geneIdtype using method if possible')
	x <- 'blabla'
	attr(x, 'annotation') <- 'toto'
	checkIdentical(fun(x), AnnotationIdentifier('toto'), "extract annotation if possible")
	x <- '10'
	attr(x, 'annotation') <- 'toto'
	checkIdentical(fun(x), EntrezIdentifier('toto'), "extract annotation if possible")
	# unigene specials
	checkIdentical(fun('Hs.1'), UnigeneIdentifier('org.Hs.eg.db'), "Extract organism from Unigene ids, if no annotation")
	checkIdentical(fun('Hs.1', annotation='toto'), UnigeneIdentifier('toto'), "Use argument annotation with Unigene ids, if present")
	# errors
	checkIdentical(fun('toto'), NullIdentifier(), "No error but return NullIdentifier() if not found")
	checkIdentical(fun('toto', annotation='aaa'), AnnotationIdentifier('aaa'), "No error but return AnnotationIdentifier() if not found but annotation is provided")
	checkException(fun('toto', error=TRUE), "error if not found and error=TRUE")
	checkIdentical(fun('toto', error='ok'), 'ok', "use value of argument error if not found")
	
}