\name{.GSE3649_pdata}
\alias{.GSE3649_pdata}
\title{Processing of Phenotypic Annotation for Dataset GSE3649}
\source{
  \url{http://genome-www.stanford.edu/normalblood/}
}
\usage{
  .GSE3649_pdata(eset, file = "GSE3649_pdata.csv")
}
\arguments{
  \item{eset}{ExpressionSet object from GSE3649 [GPL3121],
  with Nepal samples filtered out.}

  \item{file}{Phenotypic annotation file (.csv) as
  downloaded from
  \url{http://genome-www.stanford.edu/normalblood/}.}
}
\value{
  a data.frame
}
\description{
  Function used to create the additional phenotypic
  annotation data \code{\link{GSE3649_fdata}} for the
  \code{GSE3649} expression data.
}
\seealso{
  Other Whitney: \code{\link{.GSE3649_convert}}
}
\keyword{internal}

