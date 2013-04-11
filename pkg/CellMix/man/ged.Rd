\docType{methods}
\name{ged}
\alias{ged}
\alias{ged,ANY,ANY,function-method}
\alias{ged,MatrixData,ANY,ANY-method}
\alias{ged,MatrixData,ANY,GEDStrategy-method}
\alias{ged,MatrixData,ANY,missing-method}
\alias{ged-methods}
\title{Main Interface for Gene Expression Deconvolution Methods}
\usage{
  ged(object, x, method, ...)

  \S4method{ged}{MatrixData,ANY,GEDStrategy}(object, x,
    method, data = NULL, maxIter = 1, ...,
    verbose = cellmix.getOption("verbose"), wrap = TRUE,
    dryrun = FALSE)

  \S4method{ged}{ANY,ANY,function}(object, x, method, ...,
    name = NULL)
}
\arguments{
  \item{object}{global gene expression matrix-like data
  object (e.g., \code{matrix} or \code{ExpressionSet})}

  \item{x}{input data used by the algorithm to deconvolve
  global gene expression.}

  \item{method}{specification of a gene expression
  deconvolution method. If missing, a default suitable
  algorithm is selected, based on the type and/or
  dimensions of \var{x} and \var{data}.}

  \item{...}{extra arguments to allow extension, most of
  which are passed down to the deconvolution algorithm
  itself.}

  \item{data}{optional data, typically a marker list,
  specified in a format that is supported by the factory
  function \code{\link{MarkerList}}.}

  \item{maxIter}{maximum number of iterations to perform.
  If \code{method} is missing, the value of this argument
  can influence which method is selected. See section
  \emph{Details}.}

  \item{verbose}{logical that toggles verbosity. A number
  (integer) can be passed to specify the verbosity level
  (the higher the more messages are output). Passing
  \code{verbose=Inf} toggles debug mode (all messages).
  Note that because it appears after \code{...} it must be
  fully named.}

  \item{wrap}{logical that indicates the result returned by
  the method should be wrapped into an NMF object. If
  \code{FALSE}, it is returned as is. Note that because it
  appears after \code{...} it must be fully named.}

  \item{dryrun}{logical that indicates if one should make a
  fake run, i.e., without fitting any data, but showing as
  much details about the run as possible (used or inferred
  parameters, etc..). Note that because it appears after
  \code{...} it must be fully named.}

  \item{name}{optional name for the custom algorithm.}
}
\description{
  Main Interface for Gene Expression Deconvolution Methods
}
\section{Methods}{
  \describe{

  \item{ged}{\code{signature(object = "MatrixData", x =
  "ANY", method = "GEDStrategy")}: Default \code{ged}
  method apply the auto-selection scheme for determining
  which method is suitable for the type of input data. }

  \item{ged}{\code{signature(object = "MatrixData", x =
  "ANY", method = "ANY")}: Applies a deconvolution
  algorithm registered in the \pkg{CellMix} registry. }

  \item{ged}{\code{signature(object = "MatrixData", x =
  "ANY", method = "missing")}: This method deconvolves the
  target expression matrix using the algorithm selected by
  the automatic selection strategy implemented by
  \code{\link{selectGEDMethod}}, which choose a suitable
  algorithm whose data requirements match the provided
  input data (i.e. arguments \code{x} and optionally
  \code{data}). See some more details in the
  \code{\link{selectGEDMethod}} man page. }

  \item{ged}{\code{signature(object = "ANY", x = "ANY",
  method = "function")}: Applies a custom gene expression
  algorithm. }

  }
}
\examples{
# random global expression data: 3 cell types, 20 samples, 100 features
X <- rmix(3, 100, 20, markers=5)
dim(X)
# extract signature/proportion/markers
sig <- basis(X)
prop <- coef(X)
m <- getMarkers(X)
summary(m)
#--------------------------------------------
# Automatic selection of a suitable algorithm
#--------------------------------------------
# expression data only: fastdeconf
res <- ged(X, 3)
# with markers only, non iterative: qprog
res <- ged(X, m)
# with markers only, iterative: ssKL
res <- ged(X, m, maxIter=5)
# with signatures: lsfit
res <- ged(X, sig)
# with proportions: csSAM
res <- ged(X, prop)
# with proportions, iterative: DSection
if( require.quiet(RcppOctave) ){
	res <- ged(X, prop, maxIter=5)
}
}
\seealso{
  \code{\link{selectGEDMethod}}
}
\keyword{methods}

