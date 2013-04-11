\docType{methods}
\name{mixData}
\alias{mixData}
\alias{mixData<-}
\alias{mixData,ExpressionMix-method}
\alias{mixData<-,ExpressionMix,NMFstd-method}
\alias{mixData<--methods}
\alias{mixData-methods}
\title{Extracting Mixture Data}
\usage{
  mixData(object, ...)

  mixData(object) <- value
}
\arguments{
  \item{object}{object from which to extract the mixture
  data}

  \item{...}{extra argument to allow extension. See each
  method's description for more details,.}

  \item{value}{replacement value. In the case of
  \code{\link{ExpressionMix}} objects, this will be an
  object of class \code{\linkS4class{NMFstd}}.}
}
\description{
  The S4 generic function \code{mixData} extract mixture
  data from an object, i.e. information on its underlying
  constituents/components. For example, for gene expression
  deconvolution, the \pkg{CellMix} stores such data as an
  \code{\linkS4class{NMFstd}} model object.

  The function \code{mixData<-} assign mixture data onto
  objects.
}
\section{Methods}{
  \describe{

  \item{mixData}{\code{signature(object =
  "ExpressionMix")}: Extracts mixture data from an
  \code{ExpressionMix} object.

  Because \code{\linkS4class{ExpressionMix}} inherits from
  \code{linkS4class{NMFstd}}, this methods simply returns
  the the object coerced as an \code{NMFstd} object. }

  \item{mixData<-}{\code{signature(object =
  "ExpressionMix", value = "NMFstd")}: Sets mixture data in
  an ExpressionMix object.

  Because \code{\linkS4class{ExpressionMix}} inherits from
  \code{linkS4class{NMFstd}}, this method simply copies all
  slots in \code{value} in into the corresponding slots in
  \code{object}, and then checks that the result object is
  valid. }

  }
}
\examples{
e <- ExpressionMix(rmatrix(10,5))
e
# the mixture data is stored as an NMF model
mixData(e)

# assign a random mixture data with 5 cell types
m <- rnmf(3, e)
# ensure proportions sum-up to one
coef(m) <- scoef(m)
mixData(e) <- m
}
\keyword{methods}

