\docType{class}
\name{MarkerSetCollection-class}
\alias{MarkerSetCollection-class}
\title{Classes for Marker Gene Set Collections}
\description{
  This class is an extension of the
  \code{\linkS4class{GeneSetCollection}} that ensures that
  gene/feature sets are disjoint, i.e. that each
  gene/feature appear in only one set. e.g., it is used to
  store cell-specific marker genes that are used in gene
  expression deconvolution methods (cf. \code{\link{ged}}).
}
\details{
  The main difference with
  \code{\linkS4class{GeneSetCollection}} is that the
  validity method additionally ensures that there is no
  duplicated gene across all sets in the collection.
}

