\name{selectGEDMethod}
\alias{selectGEDMethod}
\title{Automatic Selection of Gene Expression Deconvolution Algorithms}
\usage{
  selectGEDMethod(object, x = NULL, data = NULL,
    maxIter = 1L, ..., call. = FALSE, quiet. = FALSE)
}
\arguments{
  \item{call.}{logical that indicates if one should return
  the suitable call to \code{\link{ged}} (\code{TRUE}), or
  just the name of the selected method.}

  \item{quiet.}{logical that indicates if an error should
  be thrown if no algorithm able to fit the input data is
  found, or simply return \code{NULL}.

  If explicitly set to \code{FALSE}, then a note is
  displayed, showing the selected algorithm and a quick
  justification for the choice.}

  \item{object}{global gene expression matrix-like data
  object (e.g., \code{matrix} or \code{ExpressionSet})}

  \item{x}{input data used by the algorithm to deconvolve
  global gene expression.}

  \item{data}{optional data, typically a marker list,
  specified in a format that is supported by the factory
  function \code{\link{MarkerList}}.}

  \item{maxIter}{maximum number of iterations to perform.
  If \code{method} is missing, the value of this argument
  can influence which method is selected. See section
  \emph{Details}.}

  \item{...}{extra arguments to allow extension, most of
  which are passed down to the deconvolution algorithm
  itself.}
}
\description{
  Implements a simple automatic selection strategy that
  chooses a suitable deconvolution method given the input
  data.
}
\details{
  The selection aims at finding an algorithm that is able
  to perform deconvolution based on the provided input
  data. The strategy is to choose amongst the possible
  algorithms available from the \pkg{CellMix}
  \emph{built-in} registry, according to their respective
  data requirements.

  Essentially the choice of algorithms made based on the
  dimensions of the target expression data \code{object}
  and the dimensions or type of the input data in \code{x}
  and \code{data}.

  Currently, the pipeline does not attempt is made to
  choose the "best" algorithm, which would be the one that
  would return the most accurate results (proportions or
  cell-specific signatures/differences) for the given data
  setting.

  The selected algorithm is indeed chosen as to be
  \emph{applicable} to the input data. When possible,
  however, a state of the art algorithm or the most
  currently used algorithm is selected.
}
\examples{
# ged methods requirements
selectGEDMethod()
# generate mixed expression data
x <- rmix(3, 100, 20)
dim(x)
sig <- basis(x)
prop <- coef(x)
ml <- getMarkers(x)

# one need at least the number of cell types
try( selectGEDMethod(x) )
selectGEDMethod(x, 3)
# from signature basis matrix
selectGEDMethod(x, sig)
selectGEDMethod(x, sig, quiet.=FALSE)
# from cell proportion matrix
selectGEDMethod(x, prop)
# from cell proportion matrix with multiple iterations
selectGEDMethod(x, prop, maxIter=10, quiet.=FALSE)
# from cell proportion matrix with markers
selectGEDMethod(x, prop, data=ml, quiet.=FALSE)
# from marker genes
selectGEDMethod(x, ml)
# from marker genes with multiple iterations
selectGEDMethod(x, ml, maxIter=10, quiet.=FALSE)
}

