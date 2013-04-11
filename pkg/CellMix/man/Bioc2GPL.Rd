\name{gpl2bioc}
\alias{bioc2gpl}
\alias{gpl2bioc}
\title{Mapping Bioconductor Annotation packages to GEO GPL Identifiers}
\usage{
  gpl2bioc(x)

  bioc2gpl(x)
}
\arguments{
  \item{x}{character vector of GPL ids or Bioconductor
  packages to map.}
}
\value{
  a character vector
}
\description{
  \code{gpl2bioc} and \code{bioc2gpl} maps a vector of GEO
  GPL ids into Bioconductor annotation packages, and vice
  versa.
}
\examples{
gpl2bioc()
gpl2bioc("GPL96")
gpl2bioc(c("GPL96", "GPL570", '123456'))
gpl2bioc("96")
bioc2gpl()
bioc2gpl('hgu133a')
bioc2gpl(c('hgu133a', 'hgu133plus2', 'abcd'))
}

