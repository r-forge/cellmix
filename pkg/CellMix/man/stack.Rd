\name{stack.MarkerList}
\alias{stack.MarkerList}
\alias{unstack.stackedMarkerList}
\title{Converting Marker Lists into data.frames}
\usage{
  \method{stack}{MarkerList} (x, ...)

  \method{unstack}{stackedMarkerList} (x, ...)
}
\arguments{
  \item{x}{a \code{\link{MarkerList}} object for
  \code{stack}, a \code{data.frame} as returned by
  \code{stack} for \code{unstack}.}

  \item{...}{extra arguments not used}
}
\description{
  The S3 method \code{stack} allow to stack a
  \code{MarkerList} object into a \code{data.frame}, which
  may be a more convenient format to handle, depending on
  the computation one wants to perform.

  The method \code{unstack} performs the reverse operation
  of \code{stack}.
}
\examples{
stack( rMarkerList(3, 3) )
stack( rMarkerList(3, 3, names=TRUE) )
stack( rMarkerList(3, 3, values=TRUE) )
m <- rMarkerList(3, 3, values=TRUE)
sm <- stack(m)
identical(unstack(sm), m)
}

