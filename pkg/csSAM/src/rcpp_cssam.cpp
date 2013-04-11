// csSAM is an R package that performs partial gene expression deconvolution
// and estimate cell-specific gene differential expression.
// Copyright (C) 2011-2012  Shai Shen-Orr, Rob Tibshirani, Narasimhan Balasubramanian, David Wang
// Copyright (C) 2012 Renaud Gaujoux (implementation in C++ of some of the functions)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "rcpp_cssam.h"

SEXP findSigGenes(SEXP rhat, SEXP cutp, SEXP fdr){
    using namespace Rcpp ;
    
    // rhat is ncell x ngenes
    NumericMatrix m_rhat(rhat);
    int numgene = m_rhat.ncol();
    int numcell = m_rhat.nrow();

    NumericMatrix m_cutp(cutp);
    NumericMatrix m_fdr(fdr);
    int thresholdLen = m_fdr.ncol();

    // initialise result matrix full of 1s
    NumericMatrix sigGene(numcell, numgene);
    std::fill(sigGene.begin(), sigGene.end(), 1);

    // find genes
    for(int curThresh=0; curThresh<thresholdLen; ++curThresh){
    		for (int curcell=0; curcell<numcell; ++curcell) {
    			for (int curgene=0; curgene<numgene; ++curgene) {
    				if (fabs(m_rhat(curcell, curgene)) >= fabs(m_cutp(curcell,curThresh))) {
    					sigGene(curcell, curgene) = m_fdr(curcell,curThresh);
    				}
    			}
    		}
    	}

    return(wrap(sigGene));
}
