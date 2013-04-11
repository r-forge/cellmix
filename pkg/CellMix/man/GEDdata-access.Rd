\docType{methods}
\name{GEDdata-access}
\alias{basis,GEDdata_entry-method}
\alias{coef,GEDdata_entry-method}
\alias{dims,GEDdata_entry-method}
\alias{exprs,GEDdata_entry-method}
\alias{GEDdata-access}
\alias{nbasis,GEDdata_entry-method}
\title{Accessing Data from CellMix Dataset}
\usage{
  \S4method{exprs}{GEDdata_entry}(object)

  \S4method{dims}{GEDdata_entry}(object)

  \S4method{basis}{GEDdata_entry}(object)

  \S4method{coef}{GEDdata_entry}(object)

  \S4method{nbasis}{GEDdata_entry}(x)
}
\arguments{
  \item{object}{an object of class \code{GEDdata_entry} as
  returned by \code{GEDdata_registry$get_entry}.}

  \item{x}{an object of class \code{"\link{factanal}"} or
  \code{"\link{princomp}"} or the \code{loadings} component
  of such an object.}
}
\value{
  a numeric matrix

  a matrix

  an integer
}
\description{
  Funtions to access data from the benchmark dataset
  registry entries
}
\section{Methods}{
  \describe{

  \item{basis}{\code{signature(object = "GEDdata_entry")}:
  \code{basis} access/compute the signatures and mixture
  Proportions. }

  \item{coef}{\code{signature(object = "GEDdata_entry")}:
  \code{coef} returns the mixture coefficient matrix (i.e.
  the proportions) available for the data. }

  \item{dims}{\code{signature(object = "GEDdata_entry")}:
  \code{dims} gets the dimensions of an ExpressionMix
  object. It returns a 3-length integer vector, containing
  the number of features, samples and components
  respectively. }

  \item{exprs}{\code{signature(object = "GEDdata_entry")}:
  \code{exprs} loads and extracts the expression matrix
  from a registry data entry. }

  \item{nbasis}{\code{signature(x = "GEDdata_entry")}:
  \code{nbasis} returns the number of constituents
  signatures available in the data. }

  }
}
\keyword{methods}

