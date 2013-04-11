\name{biocann_object}
\alias{biocann_object}
\title{Retrieving Bioconductor Annotation Maps}
\usage{
  biocann_object(to, from = NULL, optional = FALSE)
}
\arguments{
  \item{to}{target annotation field as a character string,
  e.g., \dQuote{ENTREZID}, \dQuote{ENSEMBL}, or an
  annotation package or db which means that one wants to
  retrieve a mapping to its corresponding primary
  identifier. If \code{from} is missing, \code{to} must be
  the name of an annotation package, i.e. ends with
  \dQuote{db}), in which case it tries loading the package
  and return the whole annotation db object; or any
  annotation package db or map object such as
  \code{AnnDbBimap}, \code{ChipDb} or \code{OrgDb} objects,
  which are returned unchanged.}

  \item{from}{source annotation package as a character
  string e.g. \code{"hgu133a.db"}.}

  \item{optional}{logical that indicates if the function
  should return \code{NULL} if the mapping cannot be found
  (\code{TRUE}), or throw an error. Note that this does not
  apply to the installation part: if a required annotation
  package is missing, an error is thrown even if
  \code{optional=TRUE}.}
}
\value{
  a \code{\link{ProbeAnnDbBimap}} if annotation is not
  missing, a \code{ProbeAnnDb} object otherwise.
}
\description{
  The function \code{biocann_object} retrieves annotation
  objects, like bimaps, from Bioconductor annotation
  packages. It is similar to the function \code{getAnnMap}
  from the \pkg{annotate} package, except that it also
  accepts annotation -- bimap -- objects, and will try to
  install missing packages if not found (see section
  Details).
}
\details{
  If an annotation package is specified as a character
  string, and is not found in the search path, and if R
  runs in interactive mode, then the user is asked whether
  one should try install the missing package. Default
  response is 'yes' so that simply hitting return will
  install the package via \code{\link{install.packages}}
  and load it. An error is thrown if this eventually fails.
}
\examples{
# db package object
biocann_object('hgu133a.db')

# bimap from hgu133a probe id to ENTREZID
biocann_object('ENTREZID', 'hgu133a.db')

# reversed bimap from UNIGENE to hgu133a probe id
biocann_object('hgu133a.db', 'UNIGENE')
# this is equivalent to using the annotation package object (no quotes),
# when the package is already loaded (=> helpful in interactive session with auto-completion)
biocann_object(hgu133a.db, 'UNIGENE')
}

