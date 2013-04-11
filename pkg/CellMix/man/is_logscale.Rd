\name{is_logscale}
\alias{is_logscale}
\alias{log_transform}
\title{Detect Log-transformed Data}
\source{
  \url{www.ncbi.org/geo}
}
\usage{
  is_logscale(x)

  log_transform(x, base = 2)
}
\arguments{
  \item{x}{a numeric data object (matrix, vector,
  ExpressionSet)}

  \item{base}{log base to use.}
}
\description{
  \code{is_logscale} tells if some numeric data is in log
  scale, e.g., normalized microarray data, using the same
  heuristic as GEO2R.

  \code{log_transform} apply a log transformation to the
  data. Negative values get assigned the value
  \code{\link{NaN}}.
}
\details{
  The data needs to be of reasonable size and variance for
  the detection heuristic to work correctly.
}
\examples{
x <- rmatrix(20, 10, dist=rnorm, mean=500)
is_logscale(x)
is_logscale(log_transform(x))
}

