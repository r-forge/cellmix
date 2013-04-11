\docType{methods}
\name{ExpressionMix}
\alias{ExpressionMix}
\alias{ExpressionMix,character-method}
\alias{ExpressionMix,ExpressionSet-method}
\alias{ExpressionMix,matrix-method}
\alias{ExpressionMix-methods}
\alias{show,ExpressionMix-method}
\title{Factory Method for ExpressionMix Objects}
\usage{
  ExpressionMix(object, ...)

  \S4method{show}{ExpressionMix}(object)

  \S4method{ExpressionMix}{ExpressionSet}(object, ...,
    composition = nmfModel(object))

  \S4method{ExpressionMix}{matrix}(object, ...,
    composition = nmfModel(object))

  \S4method{ExpressionMix}{character}(object,
    no.pure = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{object}{main input data. See possible values in
  section \emph{Methods}.}

  \item{...}{extra arguments to allow extension. See each
  method's description for more details.}

  \item{composition}{an \code{\linkS4class{NMF}} object,
  that contains the mixture data (signatures and/or
  proportions).}

  \item{no.pure}{a logical that indicates if the pure
  samples should be filtered out.}

  \item{verbose}{logical that indicates if log messages
  should be shown.}
}
\value{
  an \code{\linkS4class{ExpressionMix}} object
}
\description{
  This is a generic function that creates instances of
  \code{\linkS4class{ExpressionMix}} objects from common
  input data.
}
\section{Methods}{
  \describe{

  \item{ExpressionMix}{\code{signature(object =
  "ExpressionSet")}: Creates an ExpressionMix object from
  an ExpressionSet object. }

  \item{ExpressionMix}{\code{signature(object = "matrix")}:
  Creates an ExpressionMix object using the matrix
  \code{object} as the expression matrix. }

  \item{ExpressionMix}{\code{signature(object =
  "character")}: Creates an ExpressionMix object from the
  accession key of a registered dataset.

  The expression data is first loaded with
  \code{\link{gedData}} -- as an ExpressionSet object. The
  mixture data is then built as an \code{\linkS4class{NMF}}
  object and put together with the expression data into a
  new \code{\linkS4class{ExpressionMix}} object. }

  }
}
\keyword{methods}

