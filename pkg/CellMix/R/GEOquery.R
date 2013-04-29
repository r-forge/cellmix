# GEO query fixes
# 
# Author: Renaud Gaujoux
###############################################################################

## library(GEOquery)
## 
## #' Custom GEO Functions 
## #' 
## #' \code{getAndParseGSEMatrices} is an internal modified version of 
## #' \code{GEOquery::getAndParseGSEMatrices}, 
## #' which works even when the use is behind a proxy which turns GEO's dataset index 
## #' plain text pages into HTML pages.
## #' Moreover it correctly passes the destination directory argument (\code{destdir})
## #' to the subsequent calls, so that GPL description files are also downloaded there.
## #'
## #' @param GEO geo accession key
## #' @param destdir destination directory
## #' @param AnnotGPL logical that indicates if GPL data should be downloaded to 
## #' annotate the features. 
## #'
## #' @keywords internal
## #' @rdname GEO
## getAndParseGSEMatrices <- function (GEO, destdir, AnnotGPL) 
## {
##     GEO <- toupper(GEO)
##     message("# Fetching dataset index page for '", GEO, "'")
##     getURL <- get('getURL', envir = asNamespace('RCurl'))
##     a <- getURL(sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/", GEO))
##     if( grepl("^<HTML", a) ){
##         message("# Processing HTML result page (behind a proxy?) ... ", appendLF=FALSE)
##         b <- stringr::str_match_all(a, "((href)|(HREF))=\\s*[\"']/[^\"']+/([^/]+)[\"']")[[1]]
##         message('OK')
##     }else{
##         tmpcon <- textConnection(a, "r")
##         b <- read.table(tmpcon)
##         close(tmpcon)
##     }
##     b <- as.character(b[, ncol(b)])
##     message(sprintf("Found %d file(s)", length(b)))
##     ret <- list()
##     for (i in 1:length(b)) {
##         message(b[i])
##         destfile = file.path(destdir, b[i])
##         if (file.exists(destfile)) {
##             message(sprintf("Using locally cached version: %s", 
##                             destfile))
##         }
##         else {
##             download.file(sprintf("ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/%s/%s", 
##                             GEO, b[i]), destfile = destfile, mode = "wb", 
##                     method = getOption("download.file.method.GEOquery"))
##         }
##         ret[[b[i]]] <- parseGSEMatrix(destfile, AnnotGPL = AnnotGPL, destdir=destdir)$eset
##     }
##     return(ret)
## }
## 
## 
## # static variable for holding destdir for GEO
## GEOtempdir <- pkgmaker:::sVariable(NULL)
## 
## #' \code{parseGSEMatrix} fixes the issue that the destination directory for GPL 
## #' files is not passed down to other getGEO calls.
## #' 
## #' @param fname filename of the GSE matrix file
## #' 
## #' @rdname GEO
## parseGSEMatrix <- function (fname, AnnotGPL = FALSE, destdir = GEOtempdir()) 
## {
##     message("# Parsing matrix file ... ", appendLF=FALSE)
##     require(Biobase)
##     dat <- readLines(fname)
##     nseries <- sum(grepl("^!Series_", dat))
##     nsamples <- sum(grepl("^!Sample_", dat))
##     con <- GEOquery:::fileOpen(fname)
##     header <- suppressWarnings(read.table(con, sep = "\t", header = FALSE, 
##                     nrows = nseries))
##     tmpdat <- suppressWarnings(read.table(con, sep = "\t", header = FALSE, 
##                     nrows = nsamples))
##     tmptmp <- t(tmpdat)
##     sampledat <- rbind(data.frame(), tmptmp[-1, ])
##     colnames(sampledat) <- make.unique(sub("!Sample_", "", as.character(tmpdat[, 
##                                     1])))
##     suppressWarnings(readLines(con, 1))
##     datamat <- suppressWarnings(read.delim(con, sep = "\t", header = TRUE, 
##                     na.strings = c("NA", "null", "NULL", "Null"), comment.char = ""))
##     close(con)
##     tmprownames = datamat[, 1]
##     datamat <- as.matrix(datamat[!is.na(tmprownames), -1])
##     rownames(datamat) <- tmprownames[!is.na(tmprownames)]
##     datamat <- datamat[1:(nrow(datamat) - 1), ]
##     rownames(sampledat) <- colnames(datamat)
##     GPL = as.character(sampledat[1, grep("platform_id", colnames(sampledat), 
##                             ignore.case = TRUE)])
##     message("OK")
##     message("# Loading GPL data ... ", appendLF=FALSE)
##     gpl <- getGEO(GPL, AnnotGPL = AnnotGPL, destdir = destdir)
##     message("OK")
## 
##     message("# Creating annotated expression set ... ", appendLF=FALSE)
##     vmd <- Columns(gpl)
##     dat <- Table(gpl)
##     tmpnames = as.character(dat$ID)
##     tmpnames[is.na(tmpnames)] = "NA"
##     rownames(dat) <- tmpnames
##     dat <- dat[match(rownames(datamat), rownames(dat)), ]
##     rownames(vmd) <- make.unique(colnames(Table(gpl)))
##     colnames(dat) <- rownames(vmd)
##     fd <- new("AnnotatedDataFrame", data = dat, varMetadata = vmd)
##     if (is.null(nrow(datamat))) {
##         datamat = matrix(nrow = 0, ncol = nrow(sampledat))
##     }
##     eset <- new("ExpressionSet", phenoData = as(sampledat, "AnnotatedDataFrame"), 
##             annotation = GPL, featureData = fd, exprs = as.matrix(datamat))
##     message("OK")
##     return(list(GPL = as.character(sampledat[1, grep("platform_id", 
##                                             colnames(sampledat), ignore.case = TRUE)]), eset = eset))
## }

#' Downloading GSE Datasets from GEO
#' 
#' This function acts as a wrapper for \code{\link[GEOquery]{getGEO}}, to download the matrix 
#' files of gene expression datasets from the \code{GEO} database.
#' It fixes two issues:
#' 
#' \itemize{
#' \item if the user is behind a proxy that modifies the content of GEO's dataset 
#' index plain text pages into HTML pages, then plain \code{getGEO} breaks with 
#' an error like:
#' 
#' \code{Error in scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings, : line 1 did not have 6 elements}
#' \item plain \code{getGEO} do not honour the destination directory for 
#' GPL files.
#' }
#' 
#' Moreover, this function installs GEOquery if necessary.
#' 
#' @details
#' This function is a wrapper to \code{\link[GEOquery]{getGEO}}, which unlist the 
#' result list in the common case of a single dataset, and installs the \pkg{GEOquery}
#' package if necessary.
#' 
#' @inheritParams GEOquery::getGEO
#' @param ... extra parameters passed to \code{\link[GEOquery]{getGEO}}.
#' @param simplify logical that indicates if the result for GSEs that 
#' contain only one dataset should be simplified and returned as an \code{ExpressionSet} 
#' object (\code{TRUE}) or as a list.
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' dir.create('testGSE')
#' getGSE('GSE19830', destdir='testGSE')
#' }
#' 
getGSE <- function(GEO = NULL, filename = NULL, destdir = GEDtmp(), ..., simplify=TRUE)
{
	GEO <- toupper(GEO)
	# require GEOquery
	uq_requirePackage('GEOquery', load=TRUE, msg=str_c("for downloading datasets from GEO."), ptype='BioCsoft')
	
	# create a modified GEOquery "namespace"	
#	GEOns <- new.env(parent=asNamespace('GEOquery'))
#	getGEO <- GEOquery::getGEO
#	environment(getGEO) <- GEOns
#	GEOhook <- function(fname, f){
#		if( missing(f) ) f <- match.fun(fname)
#		environment(f) <- GEOns
#		GEOns[[fname]] <- f
#	}
#	GEOhook('getAndParseGSEMatrices', getAndParseGSEMatrices)
#	GEOhook('parseGSEMatrix', parseGSEMatrix)
#	GEOns$GEOtempdir <- GEOtempdir
#	GEOhook('parseGEO', GEOquery::parseGEO)
#	#
#
#	# call local getGEO
#	GEOtempdir(destdir)
#	on.exit(GEOtempdir(NULL))
	gse <- getGEO(GEO=GEO, filename=filename, destdir=destdir, ...)
	#
	
	# simplify result if necessary
	if( simplify && is.list(gse) && length(gse) == 1L ) gse[[1]]
	else gse
}


