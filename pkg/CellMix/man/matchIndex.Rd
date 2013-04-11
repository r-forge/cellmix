\docType{methods}
\name{matchIndex}
\alias{matchIndex}
\alias{matchIndex,ANY,ANY-method}
\alias{matchIndex,ANY,missing-method}
\alias{matchIndex,list,character-method}
\alias{matchIndex,list,ChipDb-method}
\alias{matchIndex,list,ExpressionSet-method}
\alias{matchIndex,list,matrix-method}
\alias{matchIndex,list,NMF-method}
\alias{matchIndex,list,ProbeAnnDbBimap-method}
\alias{matchIndex-methods}
\title{Creating Mapping from Marker Lists}
\usage{
  matchIndex(object, data, ...)

  \S4method{matchIndex}{list,character}(object, data,
    unlist = FALSE, no.match = FALSE, keep.names = TRUE,
    verbose = FALSE)
}
\arguments{
  \item{object}{to match against the reference data}

  \item{data}{reference data against which to match the
  object, i.e. retrieve the indexes from.}

  \item{...}{extra argument to allow extension}

  \item{unlist}{a logical that indicates if the indexes
  should be returned as a vector (\code{TRUE}) or another
  \code{MarkerList} object.}

  \item{no.match}{a logical that specifies if the result
  should keep the markers that cannot be found in the
  reference object. If \code{TRUE}, 0L values are used for
  unmatched markers. Default is \code{FALSE}.}

  \item{keep.names}{a logical that indicates if the
  markers' names should be kept or removed (\code{FALSE}).
  In this case, the names of the sets of each markers are
  used -- and repeated.}

  \item{verbose}{verbosity level, usually \code{TRUE} or
  \code{FALSE}.}
}
\description{
  \code{matchIndex} match indexes contained in an object
  against some other reference data, returning an object
  that is interpretable relatively to the reference data,
  e.g. retreive the row index of a data matrix contained in
  a marker list.

  Methods are provided for \code{MarkerList} objects to
  match against character vectors, \code{matrix},
  \code{\link[Biobase]{ExpressionSet-class}} or \code{NMF}
  objects, etc..
}
\section{Methods}{
  \describe{

  \item{matchIndex}{\code{signature(object = "list", data =
  "character")}: This is the workhorse method for
  \code{MarkerList} objects, that is eventually called by
  all other methods. It matches marker names against names
  in a character vector and returns the corresponding sets
  of indexes for each marker set (cell type). }

  \item{matchIndex}{\code{signature(object = "list", data =
  "matrix")}: Match markers against the row names of a data
  matrix. }

  \item{matchIndex}{\code{signature(object = "list", data =
  "NMF")}: Match markers against the basis matrix of an
  \code{NMF} model object.

  See \code{matchIndex,MarkerList,matrix-method}. }

  \item{matchIndex}{\code{signature(object = "list", data =
  "ExpressionSet")}: Match markers against the expression
  matrix of an \code{ExpressionSet} data object.

  See \code{matchIndex,MarkerList,matrix-method}. }

  \item{matchIndex}{\code{signature(object = "list", data =
  "ChipDb")}: Match markers against the keys of a chip
  annotation database. }

  \item{matchIndex}{\code{signature(object = "list", data =
  "ProbeAnnDbBimap")}: Match markers against the keys of a
  probe annotation database. }

  }
}
\examples{
# load IRIS markers and Abbas signatures
m <- cellMarkers('IRIS')
data(Abbas)
summary(m)
ma <- matchIndex(m, Abbas)
summary(ma)

# get as a mapping vector
head(matchIndex(m, Abbas, unlist=TRUE))
# only set names
head(matchIndex(m, Abbas, unlist=TRUE, keep.names=FALSE))

# if keeping all one get 0-value indexes
summary(matchIndex(m, Abbas, no.match=TRUE))
}
\keyword{methods}

