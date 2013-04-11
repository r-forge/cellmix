\docType{methods}
\name{dim,ExpressionMix-method}
\alias{dim,ExpressionMix-method}
\alias{dimnames<-,ExpressionMix-method}
\alias{dimnames,ExpressionMix-method}
\alias{featureNames<-,ExpressionMix-method}
\alias{sampleNames<-,ExpressionMix,ANY-method}
\alias{sampleNames<-,ExpressionMix-method}
\title{Dimensions in ExpressionMix Objects}
\usage{
  \S4method{dim}{ExpressionMix}(x)

  \S4method{dimnames}{ExpressionMix}(x)

  \S4method{dimnames}{ExpressionMix}(x) <- value

  \S4method{featureNames}{ExpressionMix}(object) <- value

  \S4method{sampleNames}{ExpressionMix}(object) <- value
}
\arguments{
  \item{x}{an \code{ExpressionMix} object}

  \item{value}{replacement value}

  \item{object}{Object, possibly derived from class
  \code{eSet}.}
}
\description{
  Similarly to \code{\linkS4class{NMF}} objects,
  \code{ExpressionMix} objects have a "third" dimension,
  the number of underlying cell types, inherited from the
  embedded \code{NMF} model.

  \code{dim} returns the dimensions of an ExpressionMix
  object. It returns a 3-length integer vector, containing
  the number of features, samples and components
  respectively.
}
\section{Methods}{
  \describe{

  \item{dimnames}{\code{signature(x = "ExpressionMix")}:
  \code{dimnames} returns the dimension names of an
  ExpressionMix object. It returns a list with 3 elements:
  the feature names, sample names and the constituent
  cell/tissue names (i.e. the basis names of the underlying
  NMF model). }

  \item{dimnames<-}{\code{signature(x = "ExpressionMix")}:
  \code{dimnames<-} Sets the dimension names of an
  ExpressionMix object. The replacement value must be a
  list containing the feature names, the samples names
  (optional), and the constituent cell/tissue names
  (optional), all as character vectors of correct
  dimensions (i.e. compatible with the number of features,
  samples and constituents). }

  \item{featureNames<-}{\code{signature(object =
  "ExpressionMix")}: Sets the feature names on both the
  \code{ExpressionSet} and \code{NMF} objects. }

  \item{sampleNames<-}{\code{signature(object =
  "ExpressionMix")}: Sets the sample names on both the
  \code{ExpressionSet} and \code{NMF} objects. }

  }
}
\seealso{
  \code{\link{featureNames}}, \code{\link{sampleNames}},
  \code{\link{basisnames}}
}
\keyword{methods}

