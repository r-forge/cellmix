# Unit tests for data sources
# 
# Author: Renaud Gaujoux
# Created: Apr 2, 2013
###############################################################################


DataSource <- CellMix:::DataSource

test.DataSource <- function(){
	
	# constants
	geo_ds <- 'GEO'
	ae_ds <- 'ArrayExpress'
	ds <- c('GEO', 'ArrayExpress')
	
	# checks
	checkIdentical( DataSource(),  ds, "No argument: list all datasources")
	checkIdentical( DataSource(datasource=NULL),  NULL, "Argument datasource=NULL: returns NULL.")
	lapply(ds, function(d){
		checkIdentical( DataSource(datasource=d), d, paste0("Argument datasource=", d,": returns '", d, "'"))
		checkIdentical( DataSource('GSE222', datasource=d), d, paste0("Argument datasource=", d," and valid key provided: returns '", d, "' (forced value)"))
		checkIdentical( DataSource('aaa', datasource=d), d, paste0("Argument datasource=", d," and invalid key provided: returns '", d, "' (forced value)"))
	})
	checkException( DataSource(datasource='aaa'),  "Argument datasource = invalid datasource: throw an error.")
	# GEO
	checkIdentical( DataSource('GSE124'),  geo_ds, "Argument x=GEO access key: returns GEO datasource.")
	# ArrayExpress
	lapply(c('E-TABM-124', 'E-MTAB-12', 'E-MEXP-123', 'E-GEOD-111'), function(aek){
		checkIdentical( DataSource(aek),  ae_ds, paste0("Argument x=ArrayExpress access key '", aek, "': returns ArayExpress datasource."))
	})
	TRUE
}