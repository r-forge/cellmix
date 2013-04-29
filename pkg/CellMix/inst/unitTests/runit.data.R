# Uni tets for GED datasets
# 
# Author: Renaud Gaujoux
# Creation: 13 Dec 2011
###############################################################################

library(pkgmaker)

DEACTIVATED_CRAN <- function(msg){
	if( isCRANcheck() ) DEACTIVATED(paste("CRAN check -", msg))
}

checkDataset <- function(key, edim, np, has.basis=FALSE, has.coef=FALSE, edim_eset=edim){
	
	msg <- function(...) paste(key,':', ..., sep='')
	present <- function(p){ if(p) 'present' else 'absent'} 
	d <- gedData(key)
	checkEquals(dims(d), edim, msg("GEDdata_entry dimensions are correct"), checkNames=FALSE)
	
	if( !isCRANcheck() ){
		e <- ExpressionMix(key)
		checkTrue( is(e, 'ExpressionMix'), msg("Data is an ExpressionMix object"))
		checkEquals(dim(e), edim_eset, msg("Dimensions of ExpressionMix object are corrects"), checkNames=FALSE)
		checkEquals(npure(e), np, msg("Number of pure samples is correct"), checkNames=FALSE)
		checkEquals(hasBasis(e), has.basis, msg("Basis matrix is ", present(has.basis), " as expected"))
		checkEquals(hasCoef(e), has.coef, msg("Coef matrix is ", present( has.coef ) , " as expected"))
		checkTrue( is.factor(e$Type), "Type (factor) is in pData")
		e
	}
}

# check if given proportions are present in pData
checkProp <- function(e, vars, key){	
	msg <- function(...) paste(key,':', ..., sep='')
	sapply(vars, function(t){
		checkTrue( is.numeric(e[[t]]), msg(t, " proportions are in pData"))			
	})
}

# check if cell/tissue types and their proportions are present in pData
checkType <- function(e, var, lev, key, mixname, doCheckProp=TRUE){
	
	msg <- function(...) paste(key,':', ..., sep='')
	type <- e[[var]]
	checkTrue( is.factor(type), msg(var, " is a factor in pData"))
	checkEquals( levels(type), lev, msg("Type levels are corrects"))
	
	if( doCheckProp ){
		p <- levels(type)
		if( !missing(mixname) ) # remove mixed samples
			p <- p[-grep(mixname, p)]
		checkProp(e, p, key)		
	}
}


test.GSE29832 <- function(){	
	# check dataset
	e <- checkDataset('GSE29832', c(54675, 15, 2), 6, TRUE, TRUE)	
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkType(e, 'Type', c("Blood", "Breast", "Mixed"), 'GSE29832', 'Mixed')	
}

test.GSE24223 <- function(){	
	# check dataset
	e <- checkDataset('GSE24223', c(54675, 179, 5), 110, FALSE, FALSE)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkTrue( is.factor(e$Type), "Type is in pData")
	checkEquals( levels(e$Type), c("CD14", "CD19", "CD4", "CD56", "CD8", "Whole blood"), "Type levels are corrects")	
	checkTrue( is.factor(e$Time), "Time is in pData")
	checkTrue( is.factor(e$Donor), "Donor is in pData")
}

test.GSE19830 <- function(){	
	# check dataset
	e <- checkDataset('GSE19830', c(31099, 42, 3), 9, TRUE, TRUE)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkType(e, 'Type', c("Brain", "Liver", "Lung", "Mixed"), 'GSE19830', 'Mixed')
}

test.GSE20300 <- function(){
	# check dataset
	e <- checkDataset('GSE20300', c(54675, 24, 5), 0, FALSE, TRUE)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkProp(e, c('Neutrophils', 'Lymphocytes', 'Monocytes', 'Eosinophils', 'Basophils'), 'GSE20300')
}

test.GSE5350 <- function(){
	# check dataset
	e <- checkDataset('GSE5350', c(54675, 120, 2), 60, TRUE, TRUE)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkType(e, 'Type', c("Brain", "Mixed", "Universal"), 'GSE5350', 'Mixed')
	checkTrue( is.factor(e$SampleID), "Sample is in pData")
	checkTrue( is.factor(e$Replicate), "Replicate is in pData")
	checkTrue( is.factor(e$Site), "Site is in pData")
	
}

test.GSE11057 <- function(){
	# check dataset
	e <- checkDataset('GSE11057', c(54675, 17, 3), 13, TRUE, FALSE)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkType(e, 'Type', c('Central memory', 'Effector memory', 'Naive', 'PBMCs'), 'GSE11057', 'PBMCs')
	checkTrue( is.factor(e$Donor), "Donor is in pData")
	checkTrue( is.factor(e$Group), "Group is in pData")	
	
}

test.GSE11058 <- function(){
	# check dataset
	e <- checkDataset('GSE11058', c(54675, 24, 4), 12, TRUE, TRUE)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkType(e, 'Type', sort(c('Jurkat', 'IM-9', 'Raji', 'THP-1', paste('Mix', toupper(letters[1:4]), sep=''))), 'GSE11058', '^Mix')
	checkTrue( is.factor(e$SampleID), "Sample is in pData")
	checkTrue( is.factor(e$Group), "Group is in pData")	
	
}

test.GSE22886.A <- function(){
	# check dataset
	e <- checkDataset('GSE22886_A', c(22283L, 114L, 11L), 114, FALSE, FALSE)	
}

test.GSE22886.B <- function(){
	# check dataset
	e <- checkDataset('GSE22886_B', c(22645L, 114L, 11L), 114, FALSE, FALSE)	
}

test.GSE33076 <- function(){
	# check dataset
	e <- checkDataset('GSE33076', c(22347L, 24L, 2L), 6, TRUE, TRUE)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkType(e, 'Type', c('ChAT', 'Cone', 'Mixed'), 'GSE33076', 'Mixed')
	checkTrue( is.factor(e$Line), "Line is in pData")
	checkTrue( is.factor(e$Replicate), "Replicate is in pData")
}

test.GSE3649 <- function(){
	# check dataset
	e <- checkDataset('GSE3649', c(36794L, 70L, 5L), 0, FALSE, TRUE)
}

test.Abbas <- function(){
	data(Abbas)
	checkEquals( dim(Abbas), c(359L, 17L), "Dimensions are correct", checkNames=FALSE)
}

test.SLE <- function(){
	data(SLE)
	e <- SLE
	checkEquals( dim(SLE), c(359L, 118L), "Dimensions are correct", checkNames=FALSE)
	checkTrue( is.factor(e$Type), "Type is in pData")
	checkTrue( is.factor(e$Group), "Group is in pData")
	checkIdentical( nlevels(e$Group), 2L, "Correct number of patient groups")
}

test.HaemAtlas <- function(){
	
	# modify check for buggy ArrayExpress package
	d <- c(46693L, 50L, 8L)
#	if( compareVersion(as.character(packageVersion('ArrayExpress')), "1.19.1") < 0 ){
#		d <- c(46692L, 50L, 8L)
#	}
	# check dataset
	e <- checkDataset('HaemAtlas', c(46693L, 50L, 8L), 50, TRUE, TRUE, edim_eset=d)
	
	DEACTIVATED_CRAN("Further checks require downloading dataset")
	checkType(e, 'Type', names(MarkerList('HaemAtlas')), 'HaemAtlas', doCheckProp=FALSE)
	
}
