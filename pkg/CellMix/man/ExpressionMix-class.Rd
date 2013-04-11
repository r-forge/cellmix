\docType{class}
\name{ExpressionMix-class}
\alias{ExpressionMix-class}
\title{Class for Gene Expression Deconvolution Benchmark Datasets}
\description{
  Class for Gene Expression Deconvolution Benchmark
  Datasets
}
\section{Slots}{
  \describe{

  \item{pure}{logical vector that indicates which samples
  are pure.}

  }
}
\section{Methods}{
  \describe{

  \item{dim}{\code{signature(x = "ExpressionMix")}:
  Similarly to \code{\linkS4class{NMF}} objects,
  \code{ExpressionMix} objects have a "third" dimension,
  the number of underlying cell types, inherited from the
  embedded \code{NMF} model.

  \code{dim} returns the dimensions of an ExpressionMix
  object. It returns a 3-length integer vector, containing
  the number of features, samples and components
  respectively. }

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

  \item{.DollarNames}{\code{signature(x =
  "ExpressionMix")}: Auto-completion for
  \code{\linkS4class{ExpressionMix}} objects }

  \item{eset}{\code{signature(object = "ExpressionMix")}:
  Extracts an ExpressionSet object from an
  \code{ExpressionMix} object.

  Because \code{\linkS4class{ExpressionMix}} inherits from
  \code{\linkS4class{ExpressionSet}}, this methods simply
  returns the object coerced as an \code{ExpressionSet}
  object. }

  \item{featureNames<-}{\code{signature(object =
  "ExpressionMix")}: Sets the feature names on both the
  \code{ExpressionSet} and \code{NMF} objects. }

  \item{initialize}{\code{signature(.Object =
  "ExpressionMix")}: Initialize method for ExpressionMix
  object }

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

  \item{sampleNames<-}{\code{signature(object =
  "ExpressionMix")}: Sets the sample names on both the
  \code{ExpressionSet} and \code{NMF} objects. }

  \item{showData}{\code{signature(object =
  "ExpressionMix")}: Shows which data is available in an
  \code{ExpressionMix} object.

  This method invisibly returns a 3-length logical that
  indicates if expression, signature and proportion data
  are available in \code{object}. }

  }
}

