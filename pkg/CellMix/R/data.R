################################################################################
# DATASETS DEFINITIONS
################################################################################

#' @include registry-data.R
NULL

#setGEDData('GSE38485')
#####################
# Dataset: GSE29832 [Gong et al. 2011]
#####################
gedData.GSE29832 <- setGEDData(key='GSE29832'
	, description = 'Pure/mixed blood and breast to test deconvolution of clinical samples'
	, annotation = 'hgu133plus2.db'
	, cite = 'Gong2011'
	, dim = c(54675L, 15L, 2L)
	, pure = function(eset) which(eset$Type != 'Mixed')
	, basis = meanBasis('Type')
	, coef = pdataCoef('Type')
#	function(object){
#		eset <- eset(object)
#		res <- rbind(Blood=eset$Blood, Breast=eset$Breast) / 100
#		colnames(res) <- sampleNames(eset)
#		res
#	}
	, pdata = function(eset){
		# extract type, proportions and ID
		m <- str_match(eset$title, "^([^ ]+)( \\(([0-9]+):[0-9]+\\))?,.* rep([0-9]+)")
		# type
		type <- m[,2]
		type[type=='Mixture'] <- 'Mixed'
		type <- factor(type)
		# proportions
		blood <- as.numeric(m[,4])
		blood[1:3] <- 100; blood[10:12] <- 0
		breast <- 100-blood
		biorep <- factor(m[,5])
		# return result as a data.frame
		data.frame(Type=type, Blood=blood, Breast=breast, Biorep=biorep)
	}
)

#####################
# Dataset: GSE24223
#####################

#' Dataset GSE24223: Cytometry Antigen Gene Expression
#' 
#' Supplemental Table 5 from \cite{Grigoryev2010} containing  average 
#' gene expression values for cytometry antigen markers, within each group of 
#' samples.
#' Groups are defined by the pair of experimental factors (Type, Time) (i.e. 
#' Cell type x Time of blood draw). 
#' 
#' @name GSE24223_fdata
#' @family GSE24223 
#' @source doi:10.1371/journal.pone.0013358.s007 (0.14 MB XLS)
#' \url{http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0013358.s007}
#' @keywords internal
NULL

#' \code{.GSE24223_fdata} is the function used to create the \code{\link{GSE24223_fdata}} data.
#' 
#' @param file source file
#' 
#' @source doi:10.1371/journal.pone.0013358.s007 (0.14 MB XLS)
#' \url{http://www.plosone.org/article/fetchSingleRepresentation.action?uri=info:doi/10.1371/journal.pone.0013358.s007}
#' @rdname GSE24223_fdata 
.GSE24223_fdata <- function(file="journal.pone.0013358.s007.xls"){
	# check file
	requireFile(file)
	
	# load xlsx to read Excel files
	library(xlsx)
	f <- xlsx::read.xlsx2(file, 1, startRow=2, colClasses=c(rep('character', 5), rep('numeric', 38-5)))
	# skip strange empty last column
	f <- f[,-ncol(f)]
	# remove ratio columns
	f <- f[, c(1:5, 5 + grep("average", colnames(f)[6:ncol(f)]))]
	# rename columns
	colnames(f)[1:5] <- c("PROBEID", "SYMBOL", "Antigen", "GENENAME", "Type")
	cn <- colnames(f)
	cn <- sub("Whole\\.[Bb]lood[.]+", "WB.", cn)
	cn <- sub("\\.(averaged.expression)?$", "", cn)
	colnames(f) <- cn
	# use PROBEID as rownames
	stopifnot( !anyDuplicated(f$PROBEID) )
	rownames(f) <- f$PROBEID
	# change classes
	f$Type <- factor(f$Type)
	f$Antigen <- factor(f$Antigen)
	# return data
	f
}

gedData.GSE24223 <- setGEDData(key='GSE24223'
		, description = 'Deconvoluting Early Post-Transplant Immunity Using Purified Cell Subsets'
		, annotation = 'hgu133plus2.db'
		, cite = 'Grigoryev2010'
		, dim = c(54675L, 179L, 5L)
		, pure = function(eset) which(eset$Type!='Whole blood')
		, pdata = function(eset){
			# extract donor ID
			pd <- pData(eset)[, grep("^characteristics_ch1", varLabels(eset))]
			id <- factor(str_match(pd$characteristics_ch1, "donor number: (.*)")[,2])
			# extract time point: -1=Pre-Tx, 0=Control Day 0, i:Week i
			time <- str_match(pd$characteristics_ch1.1, "time point: .* ([0-9]+)$")[,2]
			time[is.na(time)] <- -1
			time <- factor(time, levels=c(-1,0,1,2,4,8,12))
			# extract type
			type <- str_match(pd$characteristics_ch1.2, "cell subtype: (.*)")[,2]
			type[type == 'Whole Blood'] <- 'Whole blood'
			type <- factor(type)
			
			# return result as a data.frame
			data.frame(Donor=id, Time=time, Type=type)
		}
)

#####################
# Dataset: GSE19830 [Shen-Orr et al. (2010)]
#####################
gedData.GSE19830 <- setGEDData(key='GSE19830'
		, description = 'Pure/mixed brain, liver and lung to test statistical deconvolution [Rat]'
		, annotation = 'rat2302.db'
		, cite='Shen-Orr2010'
		, dim = c(31099L, 42L, 3L)
		, pure = function(eset) which(eset$Type != 'Mixed')
		, basis = meanBasis('Type')
		, coef = function(object){
			eset <- eset(object)
			t <- droplevels(eset$Type[object$pure(eset)])
			t(pData(eset)[, levels(t)]) / 100			
		}
		, pdata = function(eset){			
			# extract proportions
			comp <- eset$characteristics_ch1
			comp <- str_match_all(comp, "([0-9]+) ?% ([^0-9 ]+)")
			res <- t(sapply(comp, function(x) as.numeric(x[,2])))
			colnames(res) <-  comp[[1]][,3]
			# return result as a data.frame
			type <- apply(res, 1, function(x){
				w <- which(x==100)
				if( length(w) > 0 ) colnames(res)[w[1]]
				else "Mixed"
			})			
			res <- as.data.frame(res)
			res$Type <- factor(type)
			res
		}
)


#####################
# Dataset: GSE20300 [Shen-Orr et al. (2010)]
#####################
gedData.GSE20300 <- setGEDData(key='GSE20300'
		, description = 'Whole blood from stable and acute rejection pediatric kidney transplant'
		, annotation = 'hgu133plus2.db'
		, cite = 'Shen-Orr2010'
		, dim = c(54675L, 24L, 5L)
		, coef = function(object){
			eset <- eset(object)
			t(pData(eset)[, c('Neutrophils', 'Lymphocytes', 'Monocytes', 'Eosinophils', 'Basophils')]) / 100			
		}
		, pdata = function(eset){
			# extract state
			gr <- eset$characteristics_ch1
			gr <- str_match(gr, "state: (.*)")[,2]
			gr2 <- gr
			gr2[gr2!="Stable"] <- "ACR"
			# add phenotypic data (including proportions)
			GSE20300_pdata <- ldata('GSE20300_pdata')
			stopifnot( all(gr2 == GSE20300_pdata$Group) )
			stopifnot( all(sampleNames(eset) == GSE20300_pdata$CEL) )
			# return the result as a data.frame
			GSE20300_pdata$Status <- factor(gr)
			GSE20300_pdata$Type <- factor(rep('White blood cells',nrow(GSE20300_pdata)))
			GSE20300_pdata
		}
)

 
#####################
# Dataset: GSE5350 (MACQ) [HUMAN]
#####################
gedData.GSE5350 <- setGEDData(key='GSE5350'
		, url = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE5350/GSE5350-GPL570_series_matrix.txt.gz"
		, description = 'MicroArray Quality Control (MAQC) Project: Affymetrix HG-U133 Plus 2.0'
		, annotation = 'hgu133plus2.db'
		, cite = 'Shi2006'
		, dim = c(54675L, 120L, 2L)
		, pure = function(eset) which(eset$Type != 'Mixed') 
		, basis = meanBasis('Type')
		, coef = function(object){
			eset <- eset(object)
			t(pData(eset)[, c('Universal', 'Brain')]) / 100			
		}
		, filter = function(eset){
			ch <- as.character(eset$source_name_ch1)
			eset[, grep("^MAQC", ch)]
		}
		, pdata = function(eset){
			# extract state
			ch <- as.character(eset$source_name_ch1)
			m <- str_match(ch, "sample ([ABCD]).*((Universal)|(Brain)|(mixed))(.* at ([0-9]+)%)?")
			# sample code
			sample <- m[,2]
			# type
			type <- m[,3]
			type[type=='mixed'] <- 'Mixed'			
			# proportions
			p <- as.numeric(m[,8])
			p[type=='Universal'] <- 100
			p[type=='Brain'] <- 0
			# site: MAQC_AFX_1_A1
			ti <- eset$title
			s <- str_match(ti, "MAQC_[^_]+_([0-9]+)_")[,2]
			# replicate: MAQC_AFX_1_A1			
			re <- str_match(ti, "([0-9]+)$")[,2]
			# return as a data.frame
			data.frame(Site=factor(s), SampleID=factor(sample), Replicate=factor(re)
					, Type=factor(type), Universal=p, Brain=100-p)
		}
)

#####################
# Dataset: GSE5350 (MACQ) [RAT]
#####################
#gedData.GSE5350_Rat <- setGEDData(key='GSE5350_Rat'
#		, url = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE5350/GSE5350-GPL1355_series_matrix.txt.gz"
#		, description = 'MicroArray Quality Control (MAQC) Project: Affymetrix Rat Genome 230 2.0'
#		, dim = c(31099L, 72L, 0L)
#)

#####################
# Dataset: GSE11057 (Abbas)
#####################
gedData.GSE11057 <- setGEDData(key='GSE11057'
		, description = 'Memory T Cell Subsets: Central memory, Effector memory, Naive'
		, annotation = 'hgu133plus2.db'
		, cite = 'Abbas2009'
		, dim = c(54675L, 17L, 3L)
		, pure = function(eset) which(eset$Type!='PBMCs')
		, basis = meanBasis('Type')
#		, coef = function(object){
#			eset <- eset(object)
#			t(pData(eset)[, c('Central memory', 'Effector memory', 'Naive')]) / 100			
#		}
		, pdata = function(eset){
			# type
			type <- as.character(eset$characteristics_ch1)
			# group
			group <- type
			group[grep("PBMCs", group)] <- "Mixed"
			group[-grep("Mix", group)] <- "Pure"
			# sample id
			id <- str_match(eset$title, "([0-9]+)$")[,2]
			# proportions
			p <- NA_Matrix(ncol(eset), 3, dimnames=list(NULL, unique(type[group=='Pure'])))
			sapply(colnames(p), function(col){
				it <- type==col
				p[it,col] <<- 100
				p[it, colnames(p) != col] <<- 0
			})
			# return the result as a data.frame
			data.frame(Donor=id, Type=factor(type), Group=factor(group), as.data.frame(p), check.names=FALSE)			
		}
)

#####################
# Dataset: GSE11058 (Abbas)
#####################
gedData.GSE11058 <- setGEDData(key='GSE11058'
		, description = 'Immune Cell Line Mixtures: Jurkat, IM-9, Raji, THP-1'
		, cite = 'Abbas2009'
		, annotation = 'hgu133plus2.db'
		, dim = c(54675L, 24L, 4L)
		, pure = function(eset) which(eset$Group!='Mixed')
		, basis = meanBasis('Type')
		, coef = function(object){
			eset <- eset(object)
			t(pData(eset)[, c('Jurkat', 'IM-9', 'Raji', 'THP-1')]) / 100			
		}
		, pdata = function(eset){
			# type
			type <- as.character(eset$characteristics_ch1)
			# immune type
			ctype <- type
			ctype <- sub("(Raji)|(IM-9)","B",ctype)
			ctype <- sub("Jurkat","T",ctype)
			ctype <- sub("THP-1","Monocyte",ctype)
			# group
			group <- type
			group[grep("Mix", group)] <- "Mixed"
			group[-grep("Mix", group)] <- "Pure"
			# sample id
			id <- str_match(eset$title, "^([^ ]+)")[,2]
			# proportions
			load(packagePath('data/GSE11058_pdata.rda'), envir=environment())			
			stopifnot( all(rownames(GSE11058_pdata) == sampleNames(eset)) )
			# return the result as a data.frame
			data.frame(SampleID=id
					, Type=factor(type), CType=factor(ctype), LType=factor(paste(ctype, type, sep='-'))
					, Group=factor(group)
					, as.data.frame(GSE11058_pdata), check.names=FALSE)			
		}
)

#####################
# Dataset: GSE22886 (IRIS-A)
#####################
gedData.GSE22886_A <- setGEDData(key='GSE22886_A'
		, url = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE22886/GSE22886-GPL96_series_matrix.txt.gz"
		, description = 'IRIS: Resting and activated human immune cells [HG-U133A]'
		, cite = 'Abbas2005'
		, annotation = 'hgu133a.db'
		, dim = c(22283L, 114L, 11L)
		, pure = TRUE
		, pdata = function(eset){
			# origin tissue
			o <- str_match(eset$characteristics_ch1, "tissue: (.*)")[,2]
			# cell type
			type <- str_match(eset$characteristics_ch1.1, "type: (.*)")[,2]
			# treatment
			treat <- str_match(eset$characteristics_ch1.2, "agent: (.*)")[,2]
			treat[treat=='NA'] <- NA
			# time
			time <- str_match(eset$characteristics_ch1.3, "treatment: (.*)")[,2]			
			time[time=='NA'] <- NA
			data.frame(Origin=o, Type=type, Treatment=treat, Time=time)
		}
)

#####################
# Dataset: GSE22886 (IRIS-B)
#####################
gedData.GSE22886_B <- setGEDData(key='GSE22886_B'
		, url = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE22886/GSE22886-GPL97_series_matrix.txt.gz"
		, description = 'IRIS: Resting and activated human immune cells [HG-U133B]'
		, cite = 'Abbas2005'
		, annotation = 'hgu133b.db'
		, dim = c(22645L, 114L, 11L)
		, pure = TRUE
		, pdata = function(eset){
			# origin tissue
			o <- str_match(eset$characteristics_ch1, "tissue: (.*)")[,2]
			# cell type
			type <- str_match(eset$characteristics_ch1.1, "type: (.*)")[,2]
			# treatment
			treat <- str_match(eset$characteristics_ch1.2, "agent: (.*)")[,2]
			treat[treat=='NA'] <- NA
			# time
			time <- str_match(eset$characteristics_ch1.3, "treatment: (.*)")[,2]			
			time[time=='NA'] <- NA
			data.frame(Origin=o, Type=type, Treatment=treat, Time=time)
		}
)


#####################
# Dataset: GSE33076 [Retina mixtures] [Siegert et al. (2012)]
#####################
gedData.GSE33076 <- setGEDData(key='GSE33076'
		, description = 'Linearity of amplification between gene expression values and mRNA in retina cells'
		, cite = 'Siegert2012'
		, annotation = 'mogene10stprobeset.db'
		, dim = c(22347L, 24L, 2L)
		, pure = function(eset) which(eset$Type!='Mixed')
		, basis = meanBasis('Type')
		, coef = function(object){
			eset <- eset(object)
			t(pData(eset)[, c('ChAT', 'Cone')])	
		}
		, filter = function(eset){
			ch <- as.character(eset$characteristics_ch1)
			eset[, -grep("Retinal Arc", ch)]
		}
		, pdata = function(eset){
			# replicate
			br <- as.character(eset$title)
			br <- str_match(br, "biological rep([0-9]+)")[,2]
			# type
			ltype <- as.character(eset$characteristics_ch1)
			ltype <- str_match(ltype, "type: (.*)")[,2]
			# short type
			type <- ltype
			type[grepl('mixture', type)] <- 'Mixed'
			type[grepl('Starburst', type)] <- 'ChAT'
			type[grepl('Cone', type)] <- 'Cone'
			# proportions
			p <- as.character(eset$characteristics_ch1.2)
			p <- matrix(as.numeric(str_match(p, "ChAT ([0-9]+), cone ([0-9]+)")[,2:3])/200
						, ncol=2, dimnames=list(NULL, c('ChAT', 'Cone')))
			# line
			line <- str_match(as.character(eset$characteristics_ch1.3), "line: (.*)")[,2]
			line[grep("-RFP .*-GFP", line)] <- 'Mixed'			
			# return the result as a data.frame
			data.frame(Type=factor(type), LType=factor(ltype), Line=factor(line), Replicate=factor(br), as.data.frame(p), check.names=FALSE)			
		}
)

#####################
# Dataset: GSE3649 (Whole Blood) [Whitney et al. (2003)]
#####################
gedData.GSE3649 <- setGEDData(key='GSE3649'
		, description = 'Individuality and variation in gene expression patterns in human blood'
		, cite = 'Whitney2003'
		, url = "ftp://ftp.ncbi.nih.gov/pub/geo/DATA/SeriesMatrix/GSE3649/GSE3649-GPL3121_series_matrix.txt.gz"
		, dim = c(36794L, 70L, 5L)
		, coef = function(object){
			eset <- eset(object)
			t(pData(eset)[, c('Neutrophils', 'Lymphocytes', 'Monocytes', 'Eosinophils', 'Basophils')]) / 100	
		}
		, filter = function(eset){
			# filter out samples from Nepal which do not have phneotypic annotation
			eset <- eset[, -grep("pNepal", eset$title)]
			stopifnot( ncol(eset) == 70L )
			# filter out empty probes
			eset <- eset[fData(eset)$SPOT_ID != 'EMPTY', ]
			stopifnot( !anyMissing(fData(eset)$GB_LIST) )			
			eset
		}
#		, filter = function(gset){
#			# check that this is a list of 3 ExpressionSet objects
#			stopifnot( length(gset) == 3L )
#			# extract common probes
#			gb <- sapply(gset, function(x){
#				x <- as.character(fData(x)$GB_LIST)
#				x[ !is.na(x) & x != '' ]
#			}, simplify=FALSE)
#			# use GPL3121 as reference
#			id <- intersect(gb[[2]], gb[[3]])
#			id <- intersect(id, gb[[1]])
#			# order the data in the same order
#			for(i in 1:3) gset[[i]] <- gset[[i]][match(id, as.character(fData(gset[[i]])$GB_LIST)),]
#			# check result
#			stopifnot( all( sapply(gset, function(x) all( as.character(fData(x)$SPOT_ID) == as.character(fData(gset[[1]])$SPOT_ID)))) )
#			stopifnot( all( sapply(gset, function(x) all( varLabels(x) == varLabels(gset[[1]]) ))) )
#			
#			# combine the datasets: expression and phenotypic data
#			x <- exprs(gset[[1]])
#			pdata <- pData(gset[[1]])
#			for(i in 2:3){
#				x <- cbind(x, exprs(gset[[i]]))
#				pdata <- rbind(pdata, pData(gset[[i]]))
#			}			
#			# use spot id as feature names
#			fdata <- fData(gset[[1]])[, c('SPOT_ID', 'CONTROL', 'POLYMER', 'TYPE', 'GB_LIST')]
#			id <- as.character(fData(gset[[1]])$SPOT_ID)
#			rownames(fdata) <- id
#			rownames(x) <- id
#			# build ExpressionSet object 
#			eset <- ExpressionSet(x
#						, phenoData=AnnotatedDataFrame(pdata)
#						, featureData=AnnotatedDataFrame(fdata))
#			# filter out samples with no phenotypic annotation
#			
#		}
		, pdata = function(eset){			
			# extra annotations
			GSE3649_pdata <- ldata('GSE3649_pdata')			
			stopifnot( all(rownames(GSE3649_pdata) == sampleNames(eset)) )
			# double check consistency of the sample IDs
			sid <- str_match(as.character(eset$title), "s([0-9]+)")[,2]
			stopifnot( all(GSE3649_pdata$Sample == sid) )
			# return the result as a data.frame
			data.frame(Type="Whole blood", as.data.frame(GSE3649_pdata), check.names=FALSE)			
		}
		, fdata = function(eset){			
			# extra annotations
			GSE3649_fdata <- ldata('GSE3649_fdata')			
			stopifnot( all(rownames(GSE3649_fdata) == featureNames(eset)) )
			GSE3649_fdata			
		}
)

##########################################
# Dataset: E-TABM-633 [HaemAtlas]
##########################################
gedData.HaemAtlas <- setGEDData(key='E-TABM-633'
		, aliases = 'HaemAtlas'
		, description = 'The HaemAtlas: Transcription profiling of differentiated human blood cells'
		, cite = 'Watkins2009'
		, annotation = 'illuminaHumanv2.db'
		, dim = c(46693L, 50L, 8L)
		, pure = TRUE
		, pdata = function(eset){
			# cell type
			type <- eset$Characteristics..CellType.
			type <- .charmap(type, c(`Th lymphocyte`="T-CD4", `Tc lymphocyte`="T-CD8"
							, `monocyte`="Monocyte-CD14"
							,  `B lymphocyte`="B-CD19", `natural killer cell`="NK-CD56"
							, `granulocyte`="Granulocyte-CD66b"
							, `erythroblast`="Erythroblast", `megakaryocyte`="Megakaryocyte"))
			data.frame(Type=type)
		}
		, basis = meanBasis('Type')
		, coef = pureCoef('Type')
)

#####################
# Dataset: GSE5130
#####################
#gedData.GSE5130 <- setGEDData(key='GSE5130'
#		, description = 'Accurate and precise transcriptional profiles from 50 pg of total RNA or 100 flow-sorted primary lymphocytes'
#		, pure = integer()
#		, basis = function(object){
#			eset <- eset(object)
#			NA_Matrix(0, nrow(coef(object)))
#		} 
#		, coef = function(object){
#			eset <- eset(object)
#			t(pData(eset)[, c('T-cell', 'B-cells')]) / 100			
#		}
##		, pdata = function(eset){
##			# extract state
##			gr <- eset$characteristics_ch1
##			gr <- str_match(gr, "state: (.*)")[,2]
##			gr2 <- gr
##			gr2[gr2!="Stable"] <- "ACR"
##			load(packagePath('data/GSE20300_pdata.rda'), envir=environment())			
##			stopifnot( all(gr2 == GSE20300_pdata$Group) )
##			# return the result as a data.frame
##			GSE20300_pdata$Status <- factor(gr)
##			GSE20300_pdata
##		}
#)


################################################################################
# AUXILIARY DATA
################################################################################

#' \code{.makeData} remake some rda files shipped with CellMix. 
#' @rdname GEDdata-internals
.makeData <- function(){
	
	dir.create('data', showWarnings=FALSE)
	d <- c('Abbas', 'IRIS', 'SLE')
	sapply(d, function(x){
		message("Making data '", x, "' ... ", appendLF=FALSE)
		eval(parse(text=str_c(x, " <- .", x, "()")))
		message("OK")
		message("Saving data '", x, "' ... ", appendLF=FALSE)
		eval(parse(text=str_c("save(", x, ", file='data/", x, ".rda')")))
		message("OK")
		TRUE
	})

}
#' Cell Line Proportions for Dataset GSE11058
#' 
#' The proportions were extracted from the GEO page for dataset GSE11058.
#' 
#' @name GSE11058_pdata
#' @source \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11058}
#' @docType data
#' @cite Abbas2009
#' @family Abbas
NULL

#' Complete Blood Count for Dataset GSE20300
#' 
#' The data was extracted from the PDF file and reordered in the same order 
#' as the microarray dataset GSE20300, assuming that the GEO samples' order 
#' within each group (Stable or ACR) is consistent with the samples' order in 
#' the PDF file. 
#' 
#' @name GSE20300_pdata
#' @docType data
#' @source Extracted from Supplementary Table 3 of Shen-Orr et al. (2010)
#' @cite Shen-Orr2010
NULL


#' Basis Matrix for Whole Blood Expression Deconvolution
#' 
#' Basis matrix for performing gene expression deconvolution of whole blood 
#' as in Abbas et al. (2009).
#' 
#' The data was extracted from Supplementary Table 1 of Abbas et al. (2009), 
#' converted as text, and processed using the function \code{\link{.Abbas}}.
#' 
#' @name Abbas
#' @docType data
#' @source Supplementary Table 1 of Abbas et al. (2009)
#' Found at: doi:10.1371/journal.pone.0006098.s001 (0.11 MB PDF)
#' @cite Abbas2009
#' @family Abbas
NULL

#' Processing for Whole Blood Deconvolution Basis Matrix
#' 
#' Function used to create the deconvolution basis matrix \code{\link{Abbas}}, 
#' from supplementary data of Abbas et al. (2009).
#' 
#' @param file input file. The file would have been created from the conversion 
#' into text of the Supplementary Table S1 (a strange PDF file) from 
#' Abbas et al. (2009).
#' 
#' @cite Abbas2009
#' @family Abbas
#' @keywords internal
.Abbas <- function(file='TableS1.txt'){	
	requireFile(file, 'Abbas et al. (2009), Supplementary Table S1')
	l <- readLines(file)
	l <- str_trim(l)
	# header	
	h <- str_match_all(l[1], "([^ ]+( ((act)|([^ ]*Ig[^ ]+)))?)")[[1]]
	if( nrow(h)!=21 )
		stop("Error while parsing header of basis file: expected 20 columns [", nrow(h), "]")
	h <- h[,2]	
	# data
	ldata <- l[2:length(l)]
	d <- str_match_all(ldata, "([0-9.]+)")
	exp <- sapply(d, function(x){
				if( is.null(x) || is.na(x) || nrow(x) < 17 ) NA
				else as.numeric(tail(x[,2], 17))
	})
	if( is.list(exp) )
		stop("Error while parsing basis file: not all simple data lines return at least 17 matches")
	exp <- t(exp)
	colnames(exp) <- tail(h, ncol(exp))
	
	# extract the descriptive data
	sexp <- apply(exp, 1, paste, collapse=" ")
	l2 <- mapply(gsub, sexp, "", ldata, fixed=TRUE)
	l2 <- setNames(str_trim(l2), NULL)
	desc <- str_match(l2, "^([0-9]+[^ ]+_at)( ([0-9]+))?( ([A-Z0-9]+))?( (.+))?")
	fdata <- data.frame(AFFY=desc[,2], ENTREZID=ifelse(desc[,4]=='', NA, desc[,4])
						, SYMBOL=ifelse(desc[,6]=='', NA, desc[,6]), Name=desc[,8], stringsAsFactors=FALSE)
	
	# create phenodata	
	pdata <- read.delim(packagePath('data/Abbas.txt'), row.names=1L, sep="\t"
			, colClasses = c('character', 'character', 'factor', 'factor', 'character'))	
	
	# create ExpressionSet
	rownames(fdata) <- fdata$AFFY
	ExpressionSet(exp
				, phenoData=AnnotatedDataFrame(pdata)
				, featureData=AnnotatedDataFrame(fdata)
				, annotation = c('hgu133a.db', 'hgu133b.db'))
}

#' Gene Expression from Systemic Lupus Erythematosus Patients
#' 
#' Expression data for the basis probesets for the main clinical test cohort of 
#' 118 SLE patients and healthy controls used in Abbas et al. (2009).
#' 
#' @name SLE
#' @docType data
#' @source Extracted from Supplementary Table 2 of Abbas et al. (2009)
#' Found at: doi:10.1371/journal.pone.0006098.s002 (0.31 MB XLS)
#' @cite Abbas2009
#' @family Abbas
NULL

#' Processing for SLE Data
#' 
#' Function used to create the \code{\link{SLE}} data.
#' 
#' @param file input file. The file would be the XLS (actually tab separated) 
#' file provided as Supplementary Table S2 from Abbas et al. (2009).
#' @return an \code{\link{ExpressionSet}} object
#' @cite Abbas2009
#' @family Abbas
#' 
#' @keywords internal
.SLE <- function(file='TableS2.txt'){
	requireFile(file, 'Abbas et al. (2009), Supplementary Table S2')
	sle <- read.delim(file, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
	# extract groups
	group <- sle$group
	group[grep("Normal", group)] <- "Healthy"
	group[grep("lupus", group)] <- "SLE"
	# extract expression data
	exp <- t(sle[,-1])
	stopifnot( nrow(exp) == 359 )
	colnames(exp) <- 1:ncol(exp)
	# split into feature annotations and data
	pdata <- data.frame(Group=factor(group), Type='Whole blood')
	rownames(pdata) <- colnames(exp)
	# use feature data from the Abbas data
	Abbas <- ldata('Abbas')
	stopifnot( all(rownames(exp) == featureNames(Abbas)) )	
	# return as an ExpressionSet object
	ExpressionSet(exp
				, phenoData=AnnotatedDataFrame(pdata)
				, featureData = featureData(Abbas)
				, annotation = c('hgu133a.db', 'hgu133b.db'))	
}


#' Immune Response In-Silico (IRIS) Data
#' 
#' Extracted from Supplementary Table 1 from Abbas et al. (2005).
#' 
#' @name IRIS
#' @docType data
#' @source \url{http://share.gene.com/clark.iris.2004/irisSup/Sup6.xls}
#' @cite Abbas2005
#' @family IRIS
NULL

#' Processing for the Immune Response In-Silico (IRIS) Database
#' 
#' Function used to create the \code{\link{IRIS}} data 
#' (an \code{link{ExpressionSet}} object).
#' 
#' @param file input file. This should be a tab separated version of the 
#' Supplementary Table 1 from Abbas et al. (2005).
#' 
#' @cite Abbas2005 
#' @family IRIS
#' @keywords internal
.IRIS <- function(file='IRIS.txt'){	
	requireFile(file, 'IRIS [Abbas et al. (2005)], Supplementary Table S1')
	iris <- read.delim(file, sep="\t", stringsAsFactors=FALSE)
	rownames(iris) <- iris$Probe
	# split into feature annotations and data
	fdata <- iris[,1:11]
	colnames(fdata) <- c('AFFY','Chip','Genentech','Reference','Name', 'Description'
						,'Cell', 'NonImmuneMean', 'ImmuneMean',	'Ratio', 'F')
	fdata$Cell <- sub("[ ]*Cell", "", fdata$Cell)
	fdata$Cell <- factor(fdata$Cell)
	# add up to date annotation
	library(hgu133a.db)
	library(hgu133b.db)
	fdata$ENTREZID <- rep(NA, nrow(fdata))
	annot <- function(x, ann) sapply(AnnotationDbi::mget(fdata$AFFY[fdata$Chip==x], ann), '[', 1L)
	fdata$ENTREZID[fdata$Chip=='HGU133A'] <- annot('HGU133A', biocann_object('ENTREZID', 'hgu133a.db'))
	fdata$ENTREZID[fdata$Chip=='HGU133B'] <- annot('HGU133B', biocann_object('ENTREZID', 'hgu133b.db'))
	fdata$ENSEMBL[fdata$Chip=='HGU133A'] <- annot('HGU133A', biocann_object('ENSEMBL', 'hgu133a.db'))
	fdata$ENSEMBL[fdata$Chip=='HGU133B'] <- annot('HGU133B', biocann_object('ENSEMBL', 'hgu133b.db'))
	
	# build pheno data
	idata <- 12:33 
	cn <- colnames(iris)[idata]
	cn <- gsub("\\.\\.+", '.', cn)
	cn <- sub("T.Cells.Th", "Th", cn)
	cn <- sub("\\.Cells?", "", cn)
	cn <- sub("\\.$", "", cn)
	cn <- sub("\\.resting$", "", cn)
	cn <- sub("day\\.([0-9]+)$", "\\1d", cn)
	# type
	t <- sapply(cn, function(x){
		sapply(c(T="(\\.T)|(T\\.)|(Th[0-9]+)", NK='NK'
			, B="^B((\\.)|($))", Macrophage='Macrophage', Dendritic="Dendritic"
			, Monocyte='Monocyte', Neutrophil='Neutrophil')
		, function(p) grepl(str_c("(", p, ")"), x))		
	})
	type <- apply(t, 2L, function(x) rownames(t)[which(x)])
	stopifnot( is.vector(type) )
	# state
	state <- setNames(rep(NA, length(cn)), cn)
	state[grep("(CD[0-9]+)|(NK)|(Dendritic)|(\\.T)$", cn)] <- 'Resting'
	state[grep("(\\.act)|(IL[0-9]+)|(LPS)", cn)] <- 'Activated'
	state[grep("\\.([0-9]+h)|([0-9]+d)", cn)] <- 'Differentiation'
	pdata <- data.frame(Type=type, State=state)
	# expression data
	exp <- as.matrix(iris[,idata])	
	# sd data
	sd <- as.matrix(iris[,34:ncol(iris)])
	stopifnot( all(dim(exp)==dim(sd)) )
	rownames(pdata) <- colnames(exp) <- colnames(sd) <- cn
	
	# create ExpressionSet object
	new('ExpressionSet', list(exprs=exp, se.exprs=sd), phenoData= AnnotatedDataFrame(pdata)
					, featureData = AnnotatedDataFrame(fdata)
					, annotation = c('hgu133a.db', 'hgu133b.db'))
}

#' Phenotypic Annotation Data for Dataset GSE3649
#' 
#' This data provides several phenotypic information about the samples 
#' contained in the GEO dataset GSE3649, like age, gender, time of blood 
#' draw, as well as clinical data such as total white count, differential 
#' counts for neutrophils, lymphocytes, monocytes, eosinophils, and 
#' basophils, red blood cell count, platelet count, hemoglobin, hematocrit, 
#' and erythrocyte indices [mean corpuscular volume, mean corpuscular hemoglobin, 
#' mean corpuscular hemoglobin concentration, and red cell distribution 
#' width (RDW)].
#'    
#' Although mentioned in Whitney et al. (2003), these data were not 
#' initially available online. We thank Stephen Popper who made them available 
#' on our request.
#'  
#' @note Around November 2011, we experienced issues in accessing the web page. 
#' @name GSE3649_pdata
#' @source \url{http://genome-www.stanford.edu/normalblood/}
#' @docType data
#' @cite Whitney2003
#' @family Whitney
NULL

#' Processing of Phenotypic Annotation for Dataset GSE3649
#' 
#' Function used to create the additional phenotypic annotation data 
#' \code{\link{GSE3649_fdata}} for the \code{GSE3649} expression data.
#' 
#' @param eset ExpressionSet object from GSE3649 [GPL3121], with Nepal samples 
#' filtered out.
#' @param file Phenotypic annotation file (.csv) as downloaded from 
#' \url{http://genome-www.stanford.edu/normalblood/}.
#' @return a data.frame
#' 
#' @family Whitney
#' @source \url{http://genome-www.stanford.edu/normalblood/}
#' @keywords internal 
.GSE3649_pdata <- function(eset, file='GSE3649_pdata.csv'){
	
	requireFile(file, 'GSE3649 [Whitney et al. (2003), http://genome-www.stanford.edu/normalblood/')
	# read csv file in
	pdata <- read.csv(file, stringsAsFactors=FALSE)
	# remove duplicated columns (WBC)
	pdata <- pdata[, -grep("WBC.1", colnames(pdata))]
	# set correct class
	pdata$M.F <- factor(pdata$M.F)
	# rename some columns
	cn <- colnames(pdata)	
	ct <- c(M.F="Gender", NEUT='Neutrophils', LYM='Lymphocytes', MONO='Monocytes', EOS='Eosinophils', BASO='Basophils')
	for(i in seq_along(ct) ) cn <- str_replace_all(cn, names(ct)[i], ct[i])
	cn <- gsub(".", " ", cn, fixed=TRUE)
	colnames(pdata) <- str_trim(cn)		
	
	# add GEO sample field	
	sid <- str_match(as.character(eset$title), "s([0-9]+)")[,2]
	pid <- str_match(as.character(eset$title), "p([0-9]+)")[,2]
	stopifnot( !anyMissing(sid) && !anyMissing(pid) )
	stopifnot( !anyDuplicated(sid) && !anyDuplicated(pdata$Sample) )
	stopifnot( length(setdiff(sid, pdata$Sample)) == 0 )	
	# put in the same order as in the ExpressionSet
	pdata <- pdata[match(sid, pdata$Sample),]
	# add Patient and Array IDs
	pdata <- cbind(data.frame(CEL=sampleNames(eset), Patient=pid, stringsAsFactors=FALSE), pdata)
	rownames(pdata) <- sampleNames(eset)
	# return data.frame
	pdata
}

#' Feature Annotation Data for Dataset GSE3649
#' 
#' Extra feature annotation for the GEO dataset GSE3649.
#' 
#' We extracted and save in a file the GenBank accession numbers 
#' from the feature annotation data of the original ExpressionSet 
#' object -- created with \code{\link[GEOquery]{getGEO}}, which does not contain many 
#' annotation.
#' 
#' The content of the file (a list of GenBank IDs, one per line) was passed to 
#' IDconverter (see Reference) to obtain further annotation (ENTREZID, UNIGENE, ...).
#' We used this online service because no annotation were available for these 
#' accession numbers in the organism annotation package \code{org.Hs.eb}.  
#'   
#' @name GSE3649_fdata
#' @source \url{http://idconverter.bioinfo.cnio.es/}
#' @docType data
#' @cite Whitney2003
#' @family Whitney
NULL

#' Processing of Feature Annotation for Dataset GSE3649
#' 
#' Function used to create the additional feature annotation data 
#' \code{\link{GSE3649_fdata}} for the \code{GSE3649} expression data. 
#' 
#' The function was used twice:
#' \itemize{
#' \item once to extract and save in a file the GenBank accession numbers 
#' from the feature annotation data of the original ExpressionSet 
#' object -- created with \code{\link[GEOquery]{getGEO}}.
#' The content of the file (a list of GenBank IDs, one per line) is passed to 
#' IDconverter (see Reference) to obtain further annotation (ENTREZID, UNIGENE, ...).
#' We used this online service because no annotation were available for these 
#' accession numbers in the organism annotation package \code{org.Hs.eb}.  
#' \item once to load the annotation back from IDconverter result file (a tab 
#' delimited text file).
#' }
#' 
#' @param eset ExpressionSet object from GSE3649 [GPL3121]
#' @param file source file
#' @param mode usage mode: 'out' for writing out GenBank ids, 'in' for reading in 
#' the result file from IDconverter.
#' @return a data.frame
#' 
#' @family Whitney
#' @source \url{http://idconverter.bioinfo.cnio.es/}
#' @keywords internal 
.GSE3649_convert <- function(eset, file, mode=c('out', 'in')){
	
	mode <- match.arg(mode)
	
	if( mode == 'out '){
		# extract GenBank IDs
		gb <- fData(eset)$GB_LIST
		gb <- gb[gb!='']
		# write to file
		write(gb, file=file)
	}else{
		# load result from IDConverter
		idc <- read.delim(file, stringsAsFactors=FALSE)
		# change some column names
		newc <- c(SYMBOL="Gene", ENSEMBL="Ensembl_Gene", UNIGENE="UniGene"
				, ENTREZID="EntrezGene", ACCNUM="GenBank.Acc.", IMAGE="Clone_ID"
				, CHR="Ensembl.Chr", CHRLOC="Start..bp.", CHRLOCEND="End..bp."
				, STRAND="Strand", X="X")
		i <- match(colnames(idc), newc, nomatch=0L)
		colnames(idc)[i!=0L] <- names(newc)[i[i!=0L]]
		idc <- idc[, i!=0L]
		# reformat a bit: class, NAs where needed
		idc$ENTREZID <- as.character(idc$ENTREZID)
		sapply(c('SYMBOL', 'UNIGENE', 'ENSEMBL', 'ACCNUM', 'IMAGE', 'CHR')
				, function(x) idc[[x]][str_trim(idc[[x]])==''] <<- NA )
		# merge data
		f <- fData(eset)
		gb <- f$GB_LIST
		f <- do.call('cbind', c(list(f[,1:2]), setNames(rep(NA, ncol(idc)), colnames(idc))))
		f[gb!='', colnames(idc)] <- idc[,]
		# return the data frame
		f[, -(1:2)]
	}
	
}

.GSE3649_fdata <- function(eset, file='GSE3649_fdata.txt'){
	requireFile(file, 'GSE3649 [Whitney et al. (2003)] - GenBank ids converted by IDConverter')
	.GSE3649_convert(eset, file, mode='in')
}

