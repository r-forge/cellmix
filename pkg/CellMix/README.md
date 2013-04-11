
## Overview

This package contains methods and utilities for performing gene expression deconvolution, i.e. estimating cell type proportions and/or cell-specific gene expression signatures from global expression data in heterogeneous samples.
Its objectives are to provide:

* a unified base framework for applying and developping deconvolution analysis;
* a user-friendly and intuitive interface;
* a ressource of relevant auxiliary data such as benchmark datasets or cell marker gene lists.

## Installation
```R
# enable Bioconductor repositories
# -> add Bioc-software and Bioc-annotation
setRepositories() 

# installation
install.packages('CellMix', repos=c('http://web.cbio.uct.ac.za/~renaud/CRAN', getOption('repos')))

# or if this fails, e.g., Mac OS X, if a binary package is not available for the OS version
install.packages('CellMix', repos=c('http://web.cbio.uct.ac.za/~renaud/CRAN', getOption('repos')), type='source')
```
