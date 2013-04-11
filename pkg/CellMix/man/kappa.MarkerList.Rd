\name{kappa.MarkerList}
\alias{kappa.MarkerList}
\title{Condition Number of a Marker List}
\usage{
  \method{kappa}{MarkerList} (z, data, n = 0, ...)
}
\arguments{
  \item{z}{a \code{\linkS4class{MarkerList}} object}

  \item{data}{reference data matrix or \code{ExpressionSet}
  object, from which is extracted the sub-matrix of
  markers.}

  \item{n}{specifies the range of total number of markers
  that compose the matrix whose condition number is
  computed: \itemize{ \item The default (\code{0}) is to
  compute using all makers in \var{z} that can be found in
  the data. \item A single positive value indicates the
  exact number of such markers to use. A warning is thrown
  if this number is greater than the number of markers
  present in the data. \item A single negative value does
  the same as a positive value, but no warning is thrown if
  it exceeds the actual number of markers in the data.
  \item \code{NULL} indicates a full range of number of
  markers, i.e. the computation is performed including for
  matrices composed of \code{length(z)} to \var{n} markers.
  \item A numeric vector completely specifies the range of
  number of markers over which to perform the computation.
  }}

  \item{...}{extra arguments passed to
  \code{\link{kappa}}.}
}
\description{
  The S3 method \code{kappa} for \code{MarkerList} objects
  returns the condition number of a data (sub-)matrix
  limited to markers in the list.
}
\details{
  This method typically is used to predict/optimise the
  deconvolution power of a set of cell-specific signatures,
  as proposed by \cite{Abbas et al. (2009)}.
}
\examples{
# random data and markers
d <- rmatrix(20, 10)
x <- rMarkerList(4, 5)
summary(x)
kappa(x, d)

# over a specifc or range of total number of markers
kappa(x, d, 6)
kappa(x, d, 1:20)

# condition number is Inf if no markers are found
# NB: throws a warning
rownames(d) <- as.character(1:nrow(d))
x <- rMarkerList(4, 5, names=TRUE)
kappa(x, d)
}
\references{
  Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z and Clark
  HF (2009). "Deconvolution of blood microarray data
  identifies cellular activation patterns in systemic lupus
  erythematosus." _PloS one_, *4*(7), pp. e6098. ISSN
  1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0006098>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/19568420>.
}
\seealso{
  \code{\link{kappa}}
}

