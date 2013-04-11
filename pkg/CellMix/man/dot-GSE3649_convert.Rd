\name{.GSE3649_convert}
\alias{.GSE3649_convert}
\title{Processing of Feature Annotation for Dataset GSE3649}
\source{
  \url{http://idconverter.bioinfo.cnio.es/}
}
\usage{
  .GSE3649_convert(eset, file, mode = c("out", "in"))
}
\arguments{
  \item{eset}{ExpressionSet object from GSE3649 [GPL3121]}

  \item{file}{source file}

  \item{mode}{usage mode: 'out' for writing out GenBank
  ids, 'in' for reading in the result file from
  IDconverter.}
}
\value{
  a data.frame
}
\description{
  Function used to create the additional feature annotation
  data \code{\link{GSE3649_fdata}} for the \code{GSE3649}
  expression data.
}
\details{
  The function was used twice: \itemize{ \item once to
  extract and save in a file the GenBank accession numbers
  from the feature annotation data of the original
  ExpressionSet object -- created with
  \code{\link[GEOquery]{getGEO}}. The content of the file
  (a list of GenBank IDs, one per line) is passed to
  IDconverter (see Reference) to obtain further annotation
  (ENTREZID, UNIGENE, ...). We used this online service
  because no annotation were available for these accession
  numbers in the organism annotation package
  \code{org.Hs.eb}. \item once to load the annotation back
  from IDconverter result file (a tab delimited text file).
  }
}
\seealso{
  Other Whitney: \code{\link{.GSE3649_pdata}}
}
\keyword{internal}

