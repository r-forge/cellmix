# SCOREM approach from  Schneider et al. (2011)
# 
# SCOREM: Statistical Consolidation of Redundant Expression Measures.  NAR 2011.
#
# The main function scoremGroups is adapted from scorem::makeGroups 
# All other functions were extracted from the scorem package.
#
# Author: Renaud Gaujoux
# Created: 13 Dec 2012
###############################################################################

#' @import genefilter
#' @import graph
scoremGroups <- function (object, data, minsd=0.1, alpha=.01, method='spearman', verbose=FALSE) {
#	cat("Please be patient, this may take a few minutes.\n")
	
	eset <- data
	eset <- exprs(eset)
	# remove features that do not vary enough 
	if ( !is_NA(minsd) )  eset <- filterGenes(eset, minsd)
#	if (is.null(ann.table))  ann.table <- makeAnnTable(annotation(eset))
	
	# subset/convert marker list to match expression data
	x <- subset(object, eset, verbose=verbose)
	x <- drop(x)
	#
	
	nc <- ncol(eset)
	r.crit <- getRcrit(alpha, nc)
	w.crit <- (1+r.crit)/2
	
#	ann.table <- subset(ann.table, ann.table$probe_id %in% rownames(exprs.eset))
#	multiples <- subset(ann.table,ann.table$Freq>1)
#	new.ps.num <- by(multiples,multiples$gene_id,function(ss) {
	if( verbose ) message("Computing SCOREM groups for ", nmark(x), " markers ... ", appendLF=FALSE)
	res <- sapply(x, function(m){
		ids <- marknames(m)
		eset.rows <- eset[ids, , drop=FALSE]
		if (nrow(eset.rows)>1) { # multiple ids
			mat <- cor(t(eset.rows), method=method)
			w   <- mean(mat,na.rm=TRUE) 
			if (w < w.crit) { # find sub-groups of markers
				sg <- getSubGraphs(mat, alpha, nc, w.crit, returnScore = TRUE)  # returns names not indices
				res <- lapply(sg, function(g){
					s <- g$ids
					if (!is.numeric(s))  s <- sort(match(s,ids))
					list(ids=ids[s], score=g$score) #ss$ps.num[s] <- paste(ss$ps.num[s],collapse="_")
				})
			}else{
				res <- list(list(ids=ids, score=w))
			}
		}else res <- list(list(ids=ids, score=0))
		res# ss$ps.num
	}, simplify=FALSE)
	if( verbose ) message("OK")
	res
#	multiples$ps.num <- unlist(new.ps.num,use.names=F)
#	singles <- subset(ann.table,ann.table$Freq==1)
#	multiples <- rbind(multiples,singles)
#	new.eg.ids <- paste(multiples$gene_id,multiples$ps.num,sep='.')
#	new.eg.ids <- paste(new.eg.ids,multiples$Freq,sep='.')  
#	new.ids <- data.frame(probe_id=multiples$probe_id, new.gene_id=new.eg.ids)
#	new.ids
}

getRcrit <-	function(a, n)  {
# a is alpha, the p value cutoff
# n is the number of observations ie samples not judges (n-2 degrees of freedom)
	t <- 0-qt(a,n-2)
	sqrt(t^2/(t^2+n-2))
}

sdOverS <- function (S, na.rm = TRUE) {
	function(x) {
		sd(x, na.rm = na.rm) > S
	}
}

filterGenes <- function (eset, s=0.1) {
	ee <- exprs(eset)
	f1 <- sdOverS(s)
	ff <- filterfun(f1)
	gf <- genefilter(ee,ff)
	if( isExpressionSet(eset) ) exprs(eset) <- ee[gf,]
	else eset <- ee[gf,]
	eset
}

getSubGraphs <- function (object, alpha, nc, w.crit, returnScore=FALSE) {
	
	r.crit <- getRcrit(alpha, nc) 
	# w.crit stays the same, r.crit increases with recursion
	
	cc <- connComp(as(new("graphAM",object>r.crit),"graphNEL"))
	sg <- NULL
	lapply(cc, function(x){
		mx <- object[x,x]
		w  <- mean(mx, na.rm=TRUE)
		n  <- length(x)
		if (n==2) {
			if (w < w.crit) {
				x <- as.list(x)
				if( returnScore ){
					x <- lapply(x, function(i) list(ids=i, score=0))
				} 
			} 
			# otherwise, leave cc intact
		} else if (n==3) { 
			if (w < w.crit) { 
				rs <- c(mx[1,2],mx[1,3],mx[2,3])
				if (max(rs)>r.crit) {             # consolidate two
					pairs  <- list(c(1,2),c(1,3),c(2,3))
					tog <- pairs[[which.max(rs)]]   # pair to consolidate
					sep <- setdiff(1:3,tog)         # single to leave out
					x <- list(x[tog],x[sep])
					if( returnScore ){
						x[[1]] <- list(ids=x[[1]], score=mean(rs[tog]))
						x[[2]] <- list(ids=x[[2]], score=0)
					}
				} else {                          # don't consolidate any
					x <- as.list(x)
					if( returnScore ){
						x <- lapply(x, function(i) list(ids=i, score=0))
					}
				}
			} 
		} else { 
			if (w < w.crit) x <- getSubGraphs(mx, alpha, nc-2, w.crit, returnScore = returnScore)
			# otherwise, leave cc intact
		}
		
		if (!is.list(x))  { 
			x <- list(x) 
			if( returnScore ){
				x <- list(list(ids=unlist(x), score=w))
			}
		}  # NOTE: different than as.list(x)!!
		
		sg <<- c(sg, x)
	})
	sg
}
