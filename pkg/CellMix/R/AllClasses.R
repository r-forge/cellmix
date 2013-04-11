# Common S4 Classes
# 
# Author: Renaud Gaujoux
# Created: Mar 12, 2013
###############################################################################

#' @import Biobase
NULL

library(Biobase)

# union class for matrix-like data
setClassUnion('MatrixData', c('matrix', 'ExpressionSet'))

# tests if an object is a matrix-like data.
isMatrixData <- function(x) is(x, 'MatrixData')