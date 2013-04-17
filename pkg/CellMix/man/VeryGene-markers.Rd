\name{VeryGene-markers}
\alias{VeryGene-markers}
\title{VeryGene - Marker List for Human Tissues}
\source{
  \url{http://www.verygene.com}
}
\description{
  Marker gene list for Human tissues created from the
  \code{\link{VeryGene}} dataset. Genes are ordered
  according to their specificity index.
}
\details{
  Note that according to the VeryGene paper, tissue
  specific gene sets in this list are \strong{not
  disjoint}, and "express selectively in approximately two
  tissues on average". This can be important when used in
  deconvolution methods that require disjoint sets of
  marker genes. Duplicates accross tissue types can be
  removed using the function \code{\link{rmDuplicated}} or
  by filtering on the specificity score, e.g., \code{ml[ ml
  >= 8]}.
}
\examples{
vg <- MarkerList('VeryGene')

# marker sets are not disjoint
hasDuplicated(vg)

# show historgram of specificity scores
hist(vg)
# filter on scores
vg2 <- vg[vg >= 8, drop=TRUE]
summary(vg2)
}

