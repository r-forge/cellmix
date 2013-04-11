\docType{data}
\name{GSE3649_fdata}
\alias{GSE3649_fdata}
\title{Feature Annotation Data for Dataset GSE3649}
\source{
  \url{http://idconverter.bioinfo.cnio.es/}
}
\description{
  Extra feature annotation for the GEO dataset GSE3649.
}
\details{
  We extracted and save in a file the GenBank accession
  numbers from the feature annotation data of the original
  ExpressionSet object -- created with
  \code{\link[GEOquery]{getGEO}}, which does not contain
  many annotation.

  The content of the file (a list of GenBank IDs, one per
  line) was passed to IDconverter (see Reference) to obtain
  further annotation (ENTREZID, UNIGENE, ...). We used this
  online service because no annotation were available for
  these accession numbers in the organism annotation
  package \code{org.Hs.eb}.
}
\references{
  Whitney A, Diehn M, Popper S, Alizadeh A, Boldrick J,
  Relman D and Brown P (2003). "Individuality and variation
  in gene expression patterns in human blood." _Proceedings
  of the National Academy of Sciences of the United States
  of America_, *100*(4), pp. 1896. <URL:
  http://www.pnas.org/content/100/4/1896.short>.
}
\seealso{
  Other Whitney: \code{\link{.GSE3649_convert}},
  \code{\link{.GSE3649_pdata}}
}
\keyword{datasets}

