\name{addNames}
\alias{addNames}
\title{Generating Names}
\usage{
  addNames(object, prefix = NULL)
}
\arguments{
  \item{object}{object on which to add names}

  \item{prefix}{prefix to use for each name}
}
\value{
  an object of the same class as \code{object}
}
\description{
  Generating Names
}
\examples{
l <- list(1,2,3,4)
addNames(l)
addNames( setNames(l, letters[1:length(l)]) )
}

