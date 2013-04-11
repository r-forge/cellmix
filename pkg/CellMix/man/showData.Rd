\docType{methods}
\name{showData}
\alias{showData}
\alias{showData,ExpressionMix-method}
\alias{showData-methods}
\title{Show Data Available in an Object}
\usage{
  showData(object, ...)

  \S4method{showData}{ExpressionMix}(object, show = TRUE)
}
\arguments{
  \item{object}{object for which on wants to see the
  available data}

  \item{...}{extra arguments for future extension}

  \item{show}{logical that indicates if something is to be
  displayed on the standard output or not.}
}
\description{
  Shows which data is available in an object. For
  \code{\linkS4class{ExpressionMix}} objects, the these are
  amongst expression, cell signature and/or cell mixing
  proportion data.
}
\section{Methods}{
  \describe{

  \item{showData}{\code{signature(object =
  "ExpressionMix")}: Shows which data is available in an
  \code{ExpressionMix} object.

  This method invisibly returns a 3-length logical that
  indicates if expression, signature and proportion data
  are available in \code{object}. }

  }
}
\keyword{methods}

