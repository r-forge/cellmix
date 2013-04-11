# Test the marker list from the registry 
# 
# Author: Renaud Gaujoux
# Creation: 16 Jan 2012
###############################################################################



#' Test that one can load the marker lists
test.load <- function(){
	
	k <- cellMarkers()
	sapply(k, function(x){
		checkTrue(is(cellMarkers(x), 'MarkerList'), str_c("cellMarkers can load key '", x, "'"))
	})
	
}