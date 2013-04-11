.LOCAL_PKG_NAME <- 'CellMix'


# mask base `data` if not loading the namespace
if( !pkgmaker::isLoadingNamespace() ){
	data <- function(...){
		owd <- getwd()
		on.exit( setwd(owd) )
		setwd(packagePath())
		utils::data(...)
	}
}

cellmix <- function(){
	
	library(RcppOctave)
	o_source(packagePath('scripts/LCBMF/randcg.m'))
	o_source(packagePath('scripts/LCBMF/lcbnmf.m'))
	do.call('load_all', packageName())	
}


#stripLatex <- function(x){
#	gsub("\\\\.\\{(.)\\}", "\\1", x)
#}

# Automatic S4 Class for Registry Entries
setClassRegistry <- function(registry, Class, ...){
	
#	setClass(Class, representation, prototype, contains=character(),
#			validity, access, where, version, sealed, package,
#			S3methods = FALSE)
	
	f <- registry$get_fields()
	slots <- sapply(f, '[[', 'type', simplify=FALSE)
	defaults <- sapply(f, '[[', 'default', simplify=FALSE)
#	defaults <- defaults[names(defaults) %in% names(slots)]	
#	print(defaults)
	args <- list(Class, representation=do.call('representation', slots))
#					, prototype=do.call('prototype', defaults)
#	if( !hasArg('validity') ){
#		.validity <-
#		sapply(f, function(x){
#			if(x$is_mandatory)
#				function(object){
#					if()
#				}
#		})
#		args$validity <- function(object){
#			
#		}
#	}
	do.call('setClass', c(args, ...))
}
 
NotImplemented <- function(msg){
	stop("Not implemented - ", msg)
}
