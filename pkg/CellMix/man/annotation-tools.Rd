\docType{methods}
\name{is.annpkg}
\alias{biocann_mapname}
\alias{biocann_pkgname}
\alias{biocann_pkgobject}
\alias{is.anndb}
\alias{is.annpkg}
\alias{revmap,character-method}
\title{Annotation Tools}
\usage{
  is.annpkg(x)

  is.anndb(x)

  biocann_mapname(annotation, map, all = FALSE)

  biocann_pkgname(x, noext = FALSE)

  biocann_pkgobject(x)

  \S4method{revmap}{character}(x, ..., simplify = TRUE)
}
\arguments{
  \item{x}{an R object, either a character string or an
  annotation object.}

  \item{annotation}{names of an annotation package, with
  \dQuote{.db} suffix or not.}

  \item{map}{name of a map, e.g., \dQuote{ENTREZID}.}

  \item{all}{logical that indicates if all possible names
  should be returned, and only the simple concatenation of
  the annotation package's name without \dQuote{.db} and
  the map name.}

  \item{noext}{logical that indicates if returned package
  names should contain the extension '.db'.}

  \item{simplify}{a logical that indicates if the reverse
  mapping should be returned as a character vector whenever
  possible (\code{TRUE}), i.e. when this mapping is 1:1 --
  which can happen even if the original mapping is not 1:1.
  Otherwise the mapping is returned as a list.}

  \item{...}{extra arguments passed to
  \code{\link{revmap,list-method}}.}
}
\description{
  The \pkg{CellMix} package contains a few utility
  functions to facilitate working with Bioconductor
  annotations, which extends or enhance functions available
  in packages such as \pkg{annotate}.

  \code{is.annpkg} tells if an object is the name of an
  annotation package.

  \code{is.anndb} tells if an object is an annotation db
  object such as \code{hgu133a.db}.

  \code{biocann_mapname} returns the name of a map in an
  annotation package.

  \code{biocann_pkgname} returns the name of an annotation
  package, formated from character strings or extracted
  from annotation objects.

  \code{biocann_pkgobject} retrieve the \code{AnnotationDb}
  object for an annotation package given as a string. The
  package does not need to be already loaded but needs to
  be installed in a library that is in the search path.
}
\section{Methods}{
  \describe{

  \item{revmap}{\code{signature(x = "character")}: The
  method \code{revmap} for character vectors revert a
  mapping provided as a named character vector. }

  }
}
\examples{
# check annotation pkg name
is.annpkg('aaa.db')
is.annpkg(c('bbb.db', 'ccc.db'))
is.annpkg(c('ddd', 'eee.db'))
# check AnnotationDb object
library(hgu133a.db)
is.anndb(hgu133a.db)
# build annotation bimap object names
biocann_mapname('hgu133a.db', 'ENTREZID')
# '.db' extension is skipped
biocann_mapname('hgu133a', 'ENTREZID')
# get all possible map names
biocann_mapname('hgu133a.db', 'ENTREZID', all=TRUE)
# annotation package name
biocann_pkgname('aa')
# extract the package name from annotation objects
biocann_pkgname(hgu133a.db)
biocann_pkgname(hgu133aENTREZID)
# get AnnotationDb object
biocann_pkgobject('hgu133a') # extension can be skipped
# the package needs to be installed
try( biocann_pkgobject('aaa') )
}
\keyword{methods}

