# Defines registry of sets of markers
# 
# Author: Renaud Gaujoux
# Creation: 04 Jan 2012
###############################################################################

#' @include markers.R
NULL

##############################################################################
# REGISTRY DEFINITION
##############################################################################
#' Internal Registry for Marker Lists
#' 
#' Registry object of class \code{\link[registry:regobj]{registry}} for storing data on 
#' marker lists.
#' The data are stored as plain lists but are wrapped into \code{\linkS4class{MarkerList}} 
#' objects when retrieved with the factory function \code{\link{MarkerList}}.    
#' 
#' @keywords internal
GEDmarkers_registry <- registry::registry(registry_class = "GEDmarkers_registry"
		, entry_class = "GEDmarkers_entry")

# access key
GEDmarkers_registry$set_field("key", type="character"
		, is_key = TRUE
		, index_FUN = match_partial_ignorecase
		, validity_FUN = checkKey)

# MarkerList object
GEDmarkers_registry$set_field("data", is_mandatory = TRUE)
# Organism
GEDmarkers_registry$set_field("organism", type="character", is_mandatory=TRUE)
# IDtype
GEDmarkers_registry$set_field("idtype", type="character", default='')
# IDtype
GEDmarkers_registry$set_field("geneIdType", type="GeneIdentifierType", default=NullIdentifier())
# Annotation package
GEDmarkers_registry$set_field("annotation", type="character", is_mandatory=TRUE)
# Description
GEDmarkers_registry$set_field("shortDescription", type="character", default = '')
# Reference
GEDmarkers_registry$set_field("cite", type="character", default = '')
# PubMedIds
GEDmarkers_registry$set_field("pubMedIds", type="character", default = '')
# URLS
GEDmarkers_registry$set_field("urls", type="character", default = '')

# register using pkgmaker registry features
GEDmarkers_registry <- setPackageRegistry('markers', GEDmarkers_registry
						, description = 'Lists of marker gene for a variety of cell types'
						, entrydesc = 'marker gene list')
			
				
#' Register CellMix MarkerLists
#' 
#' \code{setMarkerList} register a marker list in the \pkg{CellMix} registry.
#' 
#' @param key accession key.
#' @param ... extra registry fields describing the marker list, 
#' or arguments passed to \code{\link[pkgmaker]{pkgreg_remove}}.
#' 
#' @export
#' @rdname GEDMarkers-internals
setMarkerList <- function(key, ...){
	setPackageRegistryEntry('markers', key, ...)
}

#' \code{removeMarkerList} removes a marker list 
#' from the registry.
#' 
#' @export
#' @rdname GEDMarkers-internals
removeMarkerList <- function(key, ...){
	pkgreg_remove('markers', key=key, ...)
}

##############################################################################
##############################################################################

#' @S3method format GEDmarkers_entry
format.GEDmarkers_entry <- function(x, ..., cite=FALSE){
	
	m <- cellMarkers(x$key)
	c(Key = x$key, Description = x$shortDescription
			, Organism = x$organism, Types=length(m)
			, Markers = nmark(m)
			, idType = if( x$idtype != '' ) x$idtype else idtype(m)
			, Annotation = .glue_anndb(x$annotation)
			, Reference = 
				if( cite ) packageReference(x$cite, short=TRUE) 
				else x$cite)
}

#' @S3method format GEDmarkers_registry
format.GEDmarkers_registry <- function(x, ...){
	rec <- x$get_entries()
	data.frame(t(sapply(rec, base::format, ...))[,-1,drop=FALSE], check.names=FALSE)	 
}

#' @S3method print GEDmarkers_registry
print.GEDmarkers_registry <- function(x, ...) print(format(x, ...))

#' @S3method xtable GEDmarkers_registry
#' @importFrom xtable xtable 
xtable.GEDmarkers_registry <- function(x, citet=FALSE, ...){
	d <- format(x, cite=FALSE)
	d$Reference <- paste("\\cite", if(citet) 't' else '', "{", d$Reference, "}", sep='')
	xtable::xtable(d, ...)
}

setOldClass('GEDmarkers_registry')
setMethod('keys', signature(x='GEDmarkers_registry')
	, function(x, keytype){
		setNames(x$get_field_entries('key'), NULL)
	}
)

#' Loading Marker Lists from Registry
#' 
#' @inheritParams gedAlgorithm
#' @param verbose logical that toggles verbose messages.
#' @param ... arguments passed to the marker list's loading function.
#' In particular argument \code{reload=TRUE} can be used to clear the cache and force 
#' reloading/recomputing the marker list from its primary data file.
#' 
#' @return a \code{GEDmarkers_entry} object or \code{NULL} (see argument \code{error}) 
#' 
#' @export
cellMarkers <- function(key, error=TRUE, verbose=FALSE, all=FALSE, ...){
	
	regobj <- cellMarkersInfo(FALSE)
	
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	
	if( missing(key) ){
		k <- keys(regobj)
		if( !all ) k <- grep("^[^.]", k, value=TRUE)
		return(k)
	}
	
	d <- pkgreg_fetch('markers', key=key, error=error, verbose=verbose)
	# exit if not found
	if( !is.list(d) ) return(d)
	
	# convert into a loading function
	res <- 
	if( is.character(d$data) ) d$data <- .cacheData(d$data)()	
	else if( is(d$data, 'MarkerList') ) d$data
	else if( is.function(d$data) ) d$data(...)
	else stop("Invalid stored `data` slot for marker list '", key, "'.")

	# make a MarkerList if necessary
	if( !is(res, 'MarkerList') ){
		res <- MarkerList(res)
	}

	# add extra information from the database
	res <- setSlots(res, d, all=FALSE)
	geneIdType(res) <- asGeneIdentifierType(if( d$idtype == '' ) idtype(res) else d$idtype, annotation=d$annotation)
	res@setName <- mkScalar(d$key)
		
	# return result
	res
}


#' \code{cellMarkersInfo} prints information about the registered marker gene lists or 
#' return the complete marker list registry, as a \code{registry} object.
#' 
#' @param show logical that indicates if the registry object should be printed (\code{FALSE}) 
#' or only returned invisibly (\code{FALSE}).
#' 
#' @rdname cellMarkers
#' @export
#' 
#' @examples 
#' # show marker lists and properties
#' cellMarkersInfo()
#' class(cellMarkersInfo(FALSE))
#' 
cellMarkersInfo <- function(show=TRUE){
	obj <- GEDmarkers_registry
	if( show ) print(obj)
	invisible(obj)
}

#####################
# Marker Lists
#####################
.cacheData <- function(x, postfun=NULL){
	.x <- x
	.md5 <- NULL
	.mList <- NULL
	function(reload = FALSE, ...){
		# initialise marker list
		md5 <- digest(list(...))
		if( reload || is.null(.mList) || md5 != .md5 ){			
			.mList <<- 
			if( is.character(.x) ){
				e <- new.env()
				data(list=.x, envir=e)
				e <- e[[.x]]
				if( !is.null(postfun) ) postfun(e)
				else e
			}else if( is.function(.x) )
				.x(...)
			.md5 <<- md5 
		}
		.mList
	}
}

###############################################################################
# DATA DEFINITIONS
###############################################################################

#' Marker Genes for Immune Cells (IRIS)
#' 
#' @description
#' \pkg{CellMix} access key: \dQuote{IRIS}
#' 
#' List of marker probesets for immune cell types from the 
#' \code{\link{IRIS}} dataset (\cite{Abbas2005}). 
#' It contains Affymetrix probesets of markers for B cells, 
#' dendritic cells, monocytes, neutrophils, NK cells and T cells, as well as 
#' for the lymphoid and myeloid lineages, and the general category of all 
#' these immune cells.
#'
#' @details 
#' The \code{MarkerList} object is generated by the internal function 
#' \code{.mIRIS}, to which arguments can be passed via 
#' \code{MarkerList('IRIS', ...)} to build lists, 
#' using different original identifiers and filtering threshold.
#' 
#' However, it is recommended to rather use the function \code{\link{convertIDs}}
#' to convert the primary identifiers into other identifiers, using up-to-date
#' annotation packages.
#' 
#' The markers are ordered in each cell type by decreasing F value (from 
#' \code{fData(IRIS)$F}), i.e. from the more to the less cell type-specific. 
#' 
#' @param id original identifier to use.
#' @param score name of the column from the \code{\link{IRIS}} dataset 
#' to use as a score.
#' @param threshold minimum value threshold for filtering on scores.
#' 
#' @source \url{http://share.gene.com/clark.iris.2004/iris/iris.html}
#' @cite Abbas2005
#' 
#' @family IRIS
#' @aliases IRIS-markers
#' 
#' @examples 
#' MarkerList('IRIS')
#' 
.mIRIS <- function(id=c('AFFY', 'ENTREZID', 'ENSEMBL'), score='F', threshold=0){
	id <- match.arg(id)
	score <- match.arg(score)
	# load IRIS data
	IRIS <- ldata('IRIS')
	# split into a list
	f <- fData(IRIS)
	# remove unnamed markers or with no score
	f <- f[!is.na(f[[id]]) & !is.na(f[[score]]) , ]
	# filter
	f <- f[f[[score]]>=threshold, ]
	# order
	f <- f[order(f[[score]], decreasing=TRUE), ]
	# split
	s <- setNames(f[[score]], f[[id]])
	res <- split(s, f$Cell)
	res <- res[c('B', 'T', 'NK', 'Dendritic', 'Monocyte', 'Neutrophil', 'Lymphoid', 'Myeloid', 'Multiple')]
	# convert into a MarkerList object
	MarkerList(res)
}

setMarkerList(key='IRIS'
							, data=.cacheData(.mIRIS)
							, organism='Human'
							, annotation=c('hgu133a.db', 'hgu133b.db')
							, shortDescription = 'Immune Response In Silico: B, T, NK and dendritic cells, monocytes and neutrophils'
							, cite = 'Abbas2005')

#' Optimised Set of Marker Genes for Immune Cells
#' 
#' @description
#' \pkg{CellMix} access key: \dQuote{Abbas}
#' 
#' Curated list of marker probesets for immune cell types from the 
#' \code{\link{Abbas}} data (\cite{Abbas2009}).
#' It contains Affymetrix probesets of markers for B cells, 
#' dendritic cells, monocytes, neutrophils, NK cells and T cells, as well as 
#' for the lymphoid and myeloid lineages, and the general category of all 
#' these immune cells.
#' 
#' @details
#' The \code{MarkerList} object is generated by the internal function \code{.mAbbas}, 
#' to which arguments can be passed via \code{MarkerList('Abbas', ...)} 
#' to build lists, using different original identifiers and filtering 
#' threshold.
#' 
#' However, it is recommended to rather use the function \code{\link{convertIDs}}
#' to convert the primary identifiers into other identifiers, using up-to-date
#' annotation packages.
#' 
#' The markers are ordered in each cell type by decreasing sparsity value,
#' i.e. from the more to the less cell type-specific, computed from the 
#' optimised signature matrix \code{\link{Abbas}}. 
#' 
#' @param sparsity sparsity threshold (see \code{\link{sparseness}})
#' @param id original identifier to use.
#' 
#' @aliases Abbas-markers
#' @source Supplementary Table 1 of Abbas et al. (2009)
#' Found at: doi:10.1371/journal.pone.0006098.s001 (0.11 MB PDF)
#' @cite Abbas2009
#' @seealso \code{\link{IRIS-markers}}
#' @family Abbas
#'  
#' @examples 
#' MarkerList('Abbas')
#' 
.mAbbas <- function(sparsity=0.8, id=c('AFFY', 'ENTREZID', 'SYMBOL')){
	id <- match.arg(id)
	# load IRIS data
	Abbas <- ldata('Abbas')
	# filter out unnamed markers
	Abbas <- Abbas[!is.na(fData(Abbas)[[id]]), ]
	# extract max expressing cell type
	s <- apply(exprs(Abbas), 1L, sparseness)
	m <- MarkerList(Abbas, values=FALSE)
	# only keep probeset with sufficiently sparse expression profile
	m <- sapply(m, function(x) sort(s[x], decreasing=TRUE), simplify=FALSE)
	m <- m[m>sparsity]
	drop(m)
}

setMarkerList(key='Abbas'
							, data=.cacheData(.mAbbas)
							, organism='Human'
							, annotation=c('hgu133a.db', 'hgu133b.db')
							, shortDescription = 'Optimised set of immune genes for deconvolution of blood samples'
							, cite = 'Abbas2009')

					
###############################			
# TissueDistributionDB
###############################
.mTDDB <- function(X, percent){
	# split into a marker list
	m <- by(X, X$Tissue, function(x){
				# only keep probeset with low percent overlap
				x <- x[x$Percent>=percent, ]
				stopifnot( !is.factor(x$Specificity) )
				m <- setNames(as.numeric(x$Specificity), x$UNIGENE)
				m[order(-x$Percent, m)]
			})
	m <- MarkerList(sapply(m, identity))
	drop(m)
}


#########
# HUMAN
#########
#' Human Tissue Specific Genes from TissueDistributionDB
#' 
#' Dataset \code{TissueDistributionDB_HS} contains the complete original
#' data obtained from \emph{TissueDistributionDB}.
#' It was processed using the internal function \code{.TissueDistributionDB_HS}.
#' 
#' @name TissueDistributionDB_HS
#' @docType data
#' @cite Kogenaru2009 
#' @family TDDB
NULL
.TissueDistributionDB_HS <- function(file='TissueDistribution_HS.txt'){
	res <- read.delim(file, sep='\t')
	# remove last column (due to badly generated tab file from TissueDistributionDB)
	res[, -ncol(res)]
	colnames(res)[2] <- 'UNIGENE'
	res$UNIGENE <- str_trim(res$UNIGENE)
	res$UNIGENE <- sub("ClusterID[ ]+", "", res$UNIGENE)
	res <- res[, c('UNIGENE', 'Tissue', 'Percent', 'Specificity')]
	rownames(res) <- res$UNIGENE
	res
}
#' Human Tissue Specific Genes from the TissueDistributionDB Database
#' 
#' @description
#' \pkg{CellMix} access key: \dQuote{TissueDistributionDB_HS}
#' 
#' Marker list for Human tissues based on Unigene EST clusters, 
#' downloaded from the TissueDistributionDB database (\cite{Kogenaru2009}).
#' 
#' @details
#' The \code{MarkerList} object is generated by the internal function \code{.mTDDB_HS}, 
#' to which arguments can be passed via \code{MarkerList('TDDB_HS', ...)} 
#' to build lists using a different filtering threshold.
#' 
#' @param percent minimum percentage expression value required for a gene to 
#' be included as a marker.
#' 
#' @aliases TDDB_HS
#' @aliases TissueDistributionDB_HS-markers
#' @cite Kogenaru2009
#' @family TDDB
#' 
#' @examples 
#' MarkerList('TDDB_HS')
#' 
.mTDDB_HS <- function(percent=99){
	# load data
	TissueDistributionDB_HS <- ldata('TissueDistributionDB_HS')
	# split into a marker list
	.mTDDB(TissueDistributionDB_HS, percent)
}


setMarkerList(key='TDDB_HS'
							, data=.cacheData(.mTDDB_HS)
							, organism='Human'
							, annotation='org.Hs.eg.db'
							, shortDescription = 'TissueDistributionDB: UniGene EST distribution profiles using Tissue Ontology from BRENDA'
							, cite = 'Kogenaru2009')


#######
# RAT
#######
#' Rat Tissue Specific Genes from TissueDistributionDB
#' 
#' Dataset \code{TissueDistributionDB_RN} contains the complete original
#' data obtained from \emph{TissueDistributionDB}.
#' It was processed using the internal function \code{.TissueDistributionDB_RN}.
#' 
#' @name TissueDistributionDB_RN
#' @docType data
#' @cite Kogenaru2009
#' @family TDDB
NULL

.TissueDistributionDB_RN <- function(file='TissueDistribution_RN.txt'){ 
	res <- read.delim(file, sep='\t')
	# remove last column (due to badly generated tab file from TissueDistributionDB)
	res[, -ncol(res)]
	colnames(res)[2] <- 'UNIGENE'
	res$UNIGENE <- str_trim(res$UNIGENE)
	res$UNIGENE <- sub("ClusterID[ ]+", "", res$UNIGENE)
	res <- res[, c('UNIGENE', 'Tissue', 'Percent', 'Specificity')]
	rownames(res) <- res$UNIGENE
	res
}

#' Rat Tissue Specific Genes from the TissueDistributionDB Database
#' 
#' @description
#' \pkg{CellMix} access key: \dQuote{TissueDistributionDB_RN}
#' 
#' Marker list for Rat tissues based on Unigene EST clusters, 
#' downloaded from the TissueDistributionDB database (\cite{Kogenaru2009}).
#' 
#' @details
#' The \code{MarkerList} object is generated by the internal function \code{.mTDDB_RN}, 
#' to which arguments can be passed via \code{MarkerList('TDDB_RN', ...)} 
#' to build lists using a different filtering threshold.
#' 
#' @param percent minimum percentage expression value required for a gene to 
#' be included as a marker.
#' 
#' @aliases TDDB_RN
#' @aliases TissueDistributionDB_RN-markers
#' @cite Kogenaru2009
#' @family TDDB
#' 
#' @examples 
#' MarkerList('TDDB_RN')
#' 
.mTDDB_RN <- function(percent=99){
	# load data
	TissueDistributionDB_RN <- ldata('TissueDistributionDB_RN')
	# split into a marker list
	.mTDDB(TissueDistributionDB_RN, percent)
}

setMarkerList(key='TDDB_RN'
							, data=.cacheData(.mTDDB_RN)
							, organism='Rat'
							, annotation='org.Rn.eg.db'
							, shortDescription = 'TissueDistributionDB: UniGene EST distribution profiles using Tissue Ontology from BRENDA'
							, cite = 'Kogenaru2009')


######################
# Palmer
######################

#' Palmer Dataset
#' 
#' Markers for B-cells, CD4+ and CD8+ T-cells, lymphocytes 
#' and granulocytes from \cite{Palmer2006}.
#' The dataset was created using the internal function \code{.Palmer}.
#'  
#' @source \url{http://www.biomedcentral.com/content/supplementary/1471-2164-7-115-S2.xls}
#' @cite Palmer2006
#' @family Palmer
#' @name Palmer
#' @docType data
NULL
.Palmer <- function(file='Palmer-Signatures.txt'){
	d <- read.delim(file, sep="\t", stringsAsFactors=FALSE)
	colnames(d)[c(1,3)] <- c('Type', 'SYMBOL')
	d$Type <- factor(d$Type)
	d
}

#' Marker List for Immune Cells - Palmer et al. (2006)
#' 
#' @description
#' \pkg{CellMix} access key: \dQuote{Palmer}
#' 
#' Marker list created from the \code{Palmer} dataset (\cite{Palmer2006}).
#' It contains markers for B-cells, CD4+ and CD8+ T-cells, lymphocytes 
#' and granulocytes.
#' 
#' @details
#' The \code{MarkerList} object is generated by the internal function \code{.mPalmer}, 
#' to which arguments can be passed via \code{MarkerList('Palmer', ...)} 
#' to build lists using different identifiers and/or filtering threshold.
#' 
#' However, it is recommended to rather use the function \code{\link{convertIDs}}
#' to convert the primary identifiers into other identifiers, using up-to-date
#' annotation packages.
#' 
#' @param id original identifier to use.
#' @param score name of the column from the \code{\link{IRIS}} dataset 
#' to use as a score.
#' @param threhold minimum value threshold for filtering on scores.
#' 
#' @family Palmer
#' @aliases Palmer-markers
#' @keywords internal
#' @examples
#' MarkerList('Palmer')
#' 
.mPalmer <- function(id=c('SYMBOL', 'LUID'), score=c('Corr', 'Fold'), threshold=0){
	id <- match.arg(id)
	score <- match.arg(score)
	# load data
	Palmer <- ldata('Palmer')
	X <- Palmer
	if( id == 'SYMBOL' ) # remove unknown genes
		X <- X[!grepl("^IMAGE:[0-9]+$", X$SYMBOL), ]
	
	# filter
	X <- X[X[[score]]>=threshold,]
	# order
	X <- X[order(X[[score]], decreasing=TRUE),]
	# split into a marker list
	s <- setNames(X[[score]], X[[id]])
	m <- split(s, X$Type)
	#str(m)
	MarkerList(m)
}


setMarkerList(key='Palmer'
							, data=.cacheData(.mPalmer)
							, organism='Human'
							, idtype='SYMBOL'
							, annotation='org.Hs.eg.db'
							, shortDescription = 'Markers for B-cells, CD4+ and CD8+ T-cells, lymphocytes and granulocytes'
							, cite = 'Palmer2006')

################################
# HaemAtlas - Watkins et al. 2009
################################

#' HaemAtlas Dataset - Immune Cells
#' 
#' The dataset contains markers for CD4+ and CD8+ T-cells, lymphocytes, 
#' monocytes, B-cells, NK cells, and granulocytes obtained from the HaemAtlas 
#' \cite{Watkins2009}.
#' 
#' The dataset was created using the internal function \code{.HaemAtlas}.
#'  
#' The probe ids provided in the supplementary data file are somehow incorrect, 
#' and do not correspond to the Illumina ids found in the illuminaHumanv2.db 
#' annotation package.
#' "Correct" probe ids were recovered from their nucelotide sequences, via their nuIDs.
#' 
#' @source \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2680378/bin/blood-2008-06-162958_TableS5.xls}
#' @family HaemAtlas
#' @cite Watkins2009
#' @name HaemAtlas
#' @aliases HaemAtlas
#' @docType data
NULL
.HaemAtlas <- function(file='blood-2008-06-162958_TableS5.xls'){
	
	lv <- lverbose(1L)
	on.exit( lverbose(lv) )
	
	# extract data from XLS file
	library(xlsx)
	vmessage("Pre-loading HaemAtlas XLS file '", file, "' ... ", appendLF=FALSE)
	wb <- loadWorkbook(file)
	sh <- getSheets(wb)
	ct <- names(sh)
	vmessage('OK [sheets: ', str_out(ct), ']')
	# rename cell types
	names(ct) <- .charmap(ct, c(`specific-CD4`="T-CD4", `specific-CD8`="T-CD8"
							, `specific-CD14`="Monocyte-CD14"
							,  `specific-CD19`="B-CD19", `specific-CD56`="NK-CD56"
							, `specific-CD66b`="Granulocyte-CD66b"
							, `specific-EB`="Erythroblast", `specific-MK`="Megakaryocyte"))
	# extract marker all data
	vmessage("Processing HaemAtlas XLS file ... ", appendLF=FALSE)
	d <- sapply(seq_along(ct), function(i){
		type <- names(ct)[i]
		sh <- ct[i]
		vmessage('Loading sheet ', sh, ' [', type, '] ... ', appendLF=FALSE)
		d <- read.xlsx(file, sheetName=sh, stringsAsFactors=FALSE)
		d$Type <- type
		vmessage('OK')
		d
	}, simplify=FALSE)
	d <- do.call('rbind', d)
	vmessage('OK')
	colnames(d) <- c('_ILMN', 'Sequence', '_ENSEMBL', '_TRANSCRIPT', '_PEPTIDE', '_SYMBOL', 'Type')
	d$Type <- factor(d$Type)
	# put NA where necessary
	sapply(c('_ENSEMBL', '_TRANSCRIPT', '_PEPTIDE', '_SYMBOL')
		, function(x) d[[x]][d[[x]]==''] <<- NA )
	
	# get correct Illumina IDs from the nuIDs using the illuminaHumanv2.db annotation package
	library(lumi)
	# load annotation for Illumina HumanWG-6 version 2 Expression BeadChip
	library(illuminaHumanv2.db)
	np <- sum(d[['_ILMN']] %in% keys(illuminaHumanv2ENTREZID))
	message("Checking number of ILMN probe ids in annotation package: ", np, '/', nrow(d))
	# compute nuIDs
	message("Converting probe ids to nuIDs ... ", appendLF=FALSE)
	nu <- lumi::seq2id(d$Sequence)
	stopifnot( !anyMissing(nu) )
	d$nuID <- nu
	message('OK [', length(nu), ']')
	
	# Auxiliary function to mget annotation and check unicity
	.mgetall <- function(keys, map){
		res <- .lookUp(keys, map)
		stopifnot( all(sapply(res, length) == 1) )
		unlist(res)
	}
	# correct Illumina probe id
	message("Recovering ILMN probe ids from nuIDs ... ", appendLF=FALSE)
	ilmn <- .mgetall(nu, revmap(illuminaHumanv2.db::illuminaHumanv2NUID))
	stopifnot( !anyMissing(ilmn) )
	message('OK')
	message("Checking duplicated ILMN probe ids ... ", sum(duplicated(ilmn)))
	#stopifnot( !anyDuplicated(ilmn) )
	d$ILMN <- ilmn
	## add other up to date annotations: SYMBOL, ENTREZID
	message("Updating SYMBOL and ENTREZIDs ... ", appendLF=FALSE)
	d$SYMBOL <- .mgetall(d$ILMN, illuminaHumanv2.db::illuminaHumanv2SYMBOL)
	d$ENTREZID <- .mgetall(d$ILMN, illuminaHumanv2.db::illuminaHumanv2ENTREZID)
	message('OK')
	##
	#pb <- sapply(seq_along(sym), function(i) is.na(sym[i]) || sym[i] == d$Symbol[i])	
	d
}

#' HaemAtlas Marker List for Immune Cells - Watkins et al. (2009)
#' 
#' @description
#' \pkg{CellMix} access key: \dQuote{HaemAtlas}
#' 
#' Marker list created from the \code{\link{HaemAtlas}} dataset (\cite{Watkins2009}).
#' It contains markers for CD4+ and CD8+ T-cells, lymphocytes, 
#' monocytes, B-cells, NK cells, and granulocytes.
#' 
#' @details
#' The \code{MarkerList} object is generated by the internal function \code{.mHaemAtlas}, 
#' to which arguments can be passed via \code{MarkerList('HaemAtlas', ...)} 
#' to build lists using different identifiers.
#' 
#' However, it is recommended to rather use the function \code{\link{convertIDs}}
#' to convert the primary identifiers into other identifiers, using up-to-date
#' annotation packages.
#' 
#' @param id original identifier to use.
#' 
#' @family HaemAtlas
#' @aliases HaemAtlas-markers
#' @examples
#' MarkerList('HaemAtlas')
#' 
.mHaemAtlas <- function(id=c('ILMN', 'nuID', 'ENTREZID', 'SYMBOL', '_ILMN', '_ENSEMBL', '_SYMBOL')){
	id <- match.arg(id)
	# load data
	Watkins <- ldata('HaemAtlas') # the data is called HaemAtlas
	# filter out unidentified markers
	Watkins <- Watkins[!is.na(Watkins[[id]]), ]
	
	# remove duplicates (check first that ids are disjoint between cell types)
	i <- sapply(levels(Watkins$Type), function(x){
		sum(Watkins$ILMN[Watkins$Type==x] %in% Watkins$ILMN[Watkins$Type!=x]) 
	})
	# we expect no cross-cell type markers
	stopifnot( all(i==0) )
	Watkins <- Watkins[!duplicated(Watkins$ILMN),]
	
	# wrap into a marker list
	MarkerList(setNames(Watkins$Type, Watkins[[id]]))
}

setMarkerList(key='HaemAtlas'
							, data=.cacheData(.mHaemAtlas)
							, organism='Human'
							, annotation='illuminaHumanv2.db'
							, shortDescription = 'HaemAtlas markers for CD4+ and CD8+ T-cells, lymphocytes, monocytes, B-cells, NK cells, and granulocytes'
							, cite = 'Watkins2009')

###########
# TIGER
###########
#' TiGER Database
#'  
#' TiGER is a database developed by the Bioinformatics Lab at Wilmer Eye 
#' Institute of Johns Hopkins University. 
#' The database contains tissue-specific gene expression profiles or expressed 
#' sequence tag (EST) data, cis-regulatory module (CRM) data, and combinatorial 
#' gene regulation data.
#' 
#' The dataset was created using the internal function \code{.TIGER}.
#'   
#' @cite Liu2008a
#' @source
#' \url{http://bioinfo.wilmer.jhu.edu/tiger}
#' \url{http://bioinfo.wilmer.jhu.edu/tiger/download/tss_spf_rsc.txt}
#' \url{http://bioinfo.wilmer.jhu.edu/tiger/download/hs2tissue-Table.txt}
#' @family TIGER
#' @name TIGER
#' @docType data
NULL
.TIGER <- function(tissuefile='TiGER-tissue.txt', ESTfile='TiGER-EST.txt', all=FALSE){
	# read tissue files in (skip first line)
	l <- readLines(tissuefile)[-1]
	m <- str_match_all(l, "\t([^\t]+)")
	hs <- str_match(l, "(Hs.[0-9]+)")[,2]
	d <- sapply(seq_along(hs), function(i) cbind(UNIGENE=hs[i], Tissue=m[[i]][,2]) )
	tissue <- do.call('rbind', d)
	# read EST file to extract the expression values	
	l <- str_c(readLines(ESTfile), collapse="#")
	raw.est <- scan(text=l, what=character(), sep=">")[-1]
	pest <- sapply(raw.est, function(x){
		hs <- str_match(x, "\t([0-9]+)\t")[,2]
		d <- strsplit(x, "#")[[1]][-1]
		d <- str_match(d, "([0-9]+)\t([0-9]+)\t([-0-9.e]+)\t([^\t]+)\t([0-9]+)\t([-0-9.e]+)")
		d <- d[!is.na(d[,1]), -c(1,6),drop=FALSE]
		colnames(d) <- c('TissueID', 'Tags', 'Score', 'Tissue', 'Log10P')
		cbind(UNIGENE=str_c('Hs.',hs), d)
	})
	pest <- do.call('rbind', pest)
	pest <- pest[!pest[,'Log10P']=='-',]
	if( all ) return(pest)
	# filter out genes using the default thresholds
	pest <- pest[str_c(pest[,'UNIGENE'],pest[,'Tissue']) %in% str_c(tissue[,'UNIGENE'],tissue[,'Tissue']), ]
	est <- as.data.frame(pest, stringsAsFactors=FALSE)
	dummy <- sapply(c('Tags', 'Score', 'Log10P'), function(x) est[[x]] <<- as.numeric(est[[x]]))
	stopifnot( !anyMissing(est) )
	est$Tissue <- factor(est$Tissue)
	est$TissueID <- factor(est$TissueID)
	est
}

#' TiGER - Human Tissue Specific Genes
#' 
#' @description
#' \pkg{CellMix} access key: \dQuote{TIGER}
#' 
#' Marker list created from the \code{\link{TIGER}} dataset (\cite{Liu2008a}).
#' The TiGER database is based on UniGene EST distribution profiles. 
#' The markers are defined by thresholds on enrichment scores and their 
#' associated log10 p-values.
#' 
#' @details
#' The \code{MarkerList} object is generated by the internal function \code{.mTIGER}, 
#' to which arguments can be passed via \code{MarkerList('TIGER', ...)} 
#' to build lists using different filtering threshold.
#' 
#' The function \code{\link{convertIDs}} can be used to convert the primary identifiers 
#' into other identifiers (e.g. ENTREZ ids or microarray probe ids), 
#' using up-to-date mappings from Bioconductor annotation packages.
#' 
#' @param score minimum score for a gene to be included in the list
#' @param log10P minimum value for -log10(P-value) for a gene to be 
#' included in the list
#' 
#' @family TIGER
#' @aliases TIGER-markers
#' 
#' @examples
#' MarkerList('TIGER')
#' 
.mTIGER <- function(score=5, log10P=3.47853){
	# load data
	TIGER <- ldata('TIGER')
	X <- TIGER
	# filter according to criteria
	X <- X[X$Score > score & X$Log10P > log10P, ]
	# order
	X  <- X[order(X$Score, decreasing=TRUE), ]
	# remove genes that flag
	# remove duplicates (check first that ids are disjoint between cell types)
	s <- setNames(X$Score, X$UNIGENE)
	m <- split(s, X$Tissue)	
	MarkerList(m)
}

setMarkerList(key='TIGER'
			, data=.cacheData(.mTIGER)
			, organism='Human'
			, annotation='org.Hs.eg.db'
			, shortDescription = 'TiGER database: based on UniGene EST distribution profiles'
			, cite = 'Liu2008a')


############
# VeryGene
############
					
#' VeryGene Dataset
#' 
#' VeryGene is developed as a curated, web-accessible centralized database for 
#' the annotation of tissue-specific/enriched genes (\cite{Yang2011}). 
#' 
#' It currently contains entries for 3960 human genes covering 128 normal 
#' tissue/cell types compiled from the expression profiling of two large 
#' microarray data sets from Liang et al. (2006) and Su et al. (2004).
#' 
#' It brings together much-needed information on preferred tissue/subcellular 
#' localization, functional annotation, pathway, mammalian phenotype, related 
#' diseases and targeting drug associated with any of these genes as a result 
#' of data integration from multiple sources. Information can be searched 
#' through gene, tissue and disease views and search result can be downloaded easily.
#' 
#' The dataset was created with the internal functions \code{.VeryGene} and 
#' \code{MySQLtoSQLite}, to import the data from the MySQL dump available online.  
#'    
#' @cite Yang2011
#' @source
#' \url{http://www.verygene.com}
#' \url{http://www.verygene.com/database.rar}
#' @family VeryGene
#' @name VeryGene
#' @docType data
#' @import DBI
#' 
NULL
.VeryGene <- function(indir='verygene/database/', dbfile="VeryGene.sqlite", verbose=TRUE){
	
	# set verbosity level
	if( !missing(verbose) ){
		ol <- lverbose(verbose)
		on.exit( lverbose(ol) )
	}
	
	# initialize a new database to a tempfile and copy some data.frame
	# from the base package into it
	file.remove(dbfile)
	con <- DBI::dbConnect("SQLite", dbname = dbfile)
	# clean up on exit in case of error
	on.exit( DBI::dbDisconnect(con), add=TRUE)
	
	sapply(c('tissue', 'entrez', 'ent_tis'), function(t){
		f <- str_c(t, '.sql')
		vmessage("Loading file '", f, "' ... ", appendLF = FALSE)
		MySQLtoSQLite(file.path(indir, f), conn=con)		
		vmessage("OK [", nrow(DBI::dbGetQuery(con, paste("select * from", t))), " entries]")
	})	
	file.info(dbfile)
	
	# extract marker gene information
	vmessage("Extract marker gene information ... ", appendLF = FALSE)
	d <- DBI::dbGetQuery(con, "select * from tissue as t, ent_tis as et, entrez as e WHERE t.tis_id = et.tis_id AND et.ent_id=e.ent_id")
	# order
	d <- d[order(d$index_score, decreasing=TRUE), ]
	# split scores
	s <- setNames(d$index_score, d$entrez)
	m <- split(s, factor(d$tissue))	
	vmessage("OK")
	#MarkerList(m)
	m
}

#' VeryGene - Marker List for Human Tissues
#' 
#' Marker gene list for Human tissues created from the \code{\link{VeryGene}} dataset.
#' Genes are ordered according to their specificity index.
#' 
#' @name VeryGene-markers
#' @source
#' \url{http://www.verygene.com}
#' @examples
#' MarkerList('VeryGene')
#' 
NULL
setMarkerList(key='VeryGene'
							, data=.cacheData('VeryGene', MarkerList)
							, organism='Human'
							, annotation='org.Hs.eg.db'
							, shortDescription = 'VeryGene database: based on two large microarray datasets'
							, cite = 'Yang2011')

#' Load a MySQL Dump for Using with SQLite
#' 
#' Converts a MySQL dump file into SQLite compatible SQL statements.
#' 
#' @note This function was developed and used to load the VeryGene database, 
#' in order to create the \code{\link{VeryGene}} marker list. 
#' 
#' @param file MySQL dump file
#' @param conn SQLite connection object
#' @return the statements as a character vector if \code{conn=NULL} or nothing  
#' 
MySQLtoSQLite <- function(file, conn=NULL){
	# load file in
	l <- readLines(file);
	# comment out unsupported commands
	## SET
	l <- gsub("^(SET [^ ]+ *= *\".*\";)$", "--\\1", l)
	## CREATE
	l <- gsub("(.*) ENGINE *=[^;]+;$", "\\1;", l)
	## INSERT: split into multiple inserts
	i.ins <- grep("INSERT INTO", l)
	ins <- character()
	if( length(i.ins) > 0 ){
		sapply(seq_along(i.ins), function(i){			
			li <- 
			if( i<length(i.ins) ) l[i.ins[i]:(i.ins[i+1]-1)] 
			else l[i.ins[i]:length(l)]
	
			lc <- str_c(li, collapse="\n")
			m <- str_match(lc, "(INSERT INTO[ \n]+[^ ]+[ \n]+\\([^)]+\\)[ \n]+VALUES)[ \n]+(.*)")			
			newins <- str_match_all(m[,3], "(\\([^\n]+\\))[;, \n]")[[1]][,2]			
			newins <- paste(m[,2], newins, ";")
			ins <<- c(ins, newins)
		})		
		l <- l[1:(i.ins[1]-1)]
	}
	query <- c(str_c(l, collapse="\n"), ins)
	cat(query, file=str_c(file, '.txt'), sep="\n")
	# send query to the database
	if( !is.null(conn) ){
		RSQLite::dbBeginTransaction(conn)
		sapply(seq_along(query), function(i, ...){
			q <- query[i]			
			tryCatch(DBI::dbSendQuery(q, ...),
			error = function(e){
				stop("Error executing: ", q, "\n", e$msg)
			})
		}, conn=conn)		
		DBI::dbCommit(conn)		
	}else query
} 
