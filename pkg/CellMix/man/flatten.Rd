\docType{methods}
\name{flatten}
\alias{flatten}
\alias{flatten,list-method}
\alias{flatten,MarkerList-method}
\alias{flatten-methods}
\title{Flattening Marker Lists}
\usage{
  flatten(object, ...)

  \S4method{flatten}{list}(object, use.names = FALSE)

  \S4method{flatten}{MarkerList}(object, use.names = FALSE,
    size = NA)
}
\arguments{
  \item{object}{the object to be flattened}

  \item{...}{extra arguments to allow extension. See each
  method's description for more details.}

  \item{use.names}{logical that indicates if the top level
  names, i.e. the names of list (e.g., the cell types)
  should be used instead of the names of each vector
  element.}

  \item{size}{an integer that specifies the desired size of
  the result vector, which is used only if the marker are
  identified by their index. It is useful if one wants to
  use the result with data of known -- greater --
  dimension.}
}
\description{
  \code{flatten} is an S4 generic function that flattens
  objects by unfolding their inner levels, like what
  \code{\link{unlist}} does for \code{\link{list}} objects.
}
\section{Methods}{
  \describe{

  \item{flatten}{\code{signature(object = "list")}: For
  \code{\link{list}} objects this is equivalent to
  \code{\link{unlist2}}. }

  \item{flatten}{\code{signature(object = "MarkerList")}:
  For \code{\link{MarkerList}} objects the conversion is
  similar to what \code{\link{unlist}} would do, but
  argument \code{use.names} is used slightly differently.

  Secondly, the names of the main elements are append to
  the marker names (like \code{unlist} does if
  \code{use.names=TRUE}) only for marker lists that contain
  numeric values (e.g. specificity scores). In the other
  cases, they are used either as values for character
  marker lists or as duplicated names for integer marker
  lists. }

  }
}
\examples{
l <- list(1,2,3,c(4,5,6))
flatten(l)
nl <- setNames(l, letters[1:length(l)])
flatten(nl)
flatten(nl, use.names=TRUE)
ml <- rMarkerList(3, names=TRUE)
flatten(ml)
flatten(ml, use.names=TRUE)
}
\seealso{
  \code{\link{unlist}}
}
\keyword{methods}

