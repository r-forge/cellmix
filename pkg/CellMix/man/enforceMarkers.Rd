\docType{methods}
\name{enforceMarkers}
\alias{enforceMarkers}
\alias{enforceMarkers,matrix,list-method}
\alias{enforceMarkers,matrix,numeric-method}
\alias{enforceMarkers-methods}
\alias{enforceMarkers,NMF,ANY-method}
\title{Enforcing Marker Block Patterns}
\usage{
  enforceMarkers(object, markers, ...)

  \S4method{enforceMarkers}{matrix,list}(object, markers,
    ratio = NULL, value = NULL, eps = NULL, attach = FALSE)
}
\arguments{
  \item{object}{a numeric matrix or an
  \code{\link[NMF]{NMF-class}} object}

  \item{markers}{a numeric vector or the list of markers to
  enforce. This is usually an \code{MarkerList} object, but
  can be any list, as long as it contains markers that can
  matched to rows in \code{object}. If \code{markers} is a
  single numeric, then it specifies the number of markers
  to enforce in each cell-type. The markers are enforced on
  the top -- basis -- rows of \code{object}. If
  \code{markers} is a numeric vector, then its length must
  be the number of columns (basis) in \code{object}.}

  \item{ratio}{a single numeric value specifying the
  minimum expression fold change to enforce. If \code{NULL}
  or zero, then the expression values of each marker are
  forced to 0 on non-related cell types. Default is
  \code{NULL}.}

  \item{value}{a single numeric value specifying the
  expression value to enforce on the cell-type related to
  each marker. If \code{NULL} this value is left unchanged.
  Default is \code{NULL}.}

  \item{eps}{minimum value imposed for null entries.}

  \item{attach}{logical that indicates if the markers
  should be attached to the matrix data, as attibute
  \code{'markers'}.}

  \item{...}{extra arguments to allow extension, and passed
  down to the workhorse method
  \code{enforceMarkers,matrix,list}.}
}
\value{
  the input object with the marker block patterns enforced
}
\description{
  This method enforces marker block patterns on a numeric
  matrix, which typically is the basis matrix of the NMF
  model being estimated.
}
\details{
  The block patterns are defined by a list of markers whose
  expression is cell-type specific. The expression value
  \eqn{e_{ij}}{e_ij} for a given marker \eqn{i}, on its
  respective cell-type \eqn{j} is forced to a given value;
  while its expression \eqn{e_{ik}}{e_ik} on any other
  cell-type \eqn{k} is forced to 0 or to be a certain times
  lower than \eqn{e_{ij}}{e_ij}: \eqn{e_{ik} =
  \min(\frac{e_{ij}}{f} - \epsilon, e_{ik})}{e_ik =
  min(e_ij/f - EPS, e_ik)}.
}
\section{Methods}{
  \describe{

  \item{enforceMarkers}{\code{signature(object = "matrix",
  markers = "list")}: Enforce marker patterns specified as
  a list.

  This is the workhorse method that is eventually called by
  all other methods. }

  \item{enforceMarkers}{\code{signature(object = "matrix",
  markers = "numeric")}: Enforce a given number of marker
  patterns }

  \item{enforceMarkers}{\code{signature(object = "NMF",
  markers = "ANY")}: The method for an \code{NMF} object
  enforces the markers on its basis matrix and sets its
  basis names to the marker list's names, if none were
  already defined. }

  }
}
\author{
  Renaud Gaujoux
}
\keyword{methods}

