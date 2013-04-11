\name{gedProportions}
\alias{gedProportions}
\title{Estimating Cell Proportions from Known Signatures}
\usage{
  gedProportions(object, x, method = "lsfit",
    subset = NULL,
    map.method = if (is.matrix(x)) "names" else "auto",
    ..., log = NULL, lbase = 2,
    normalize = c("none", "quantiles"), verbose = FALSE)
}
\arguments{
  \item{object}{target data, specified in any format
  supported by \code{\link{ged}}.}

  \item{x}{basis signature data. It can be an
  \code{\linkS4class{ExpressionSet}} object or a
  \code{matrix}, whose columns contains cell-specific
  expression for each feature in the target data.

  If the gene identifier type from the basis matrix do not
  match the one from the target matrix, these are converted
  using \code{\link{convertIDs}}. If needed, this automatic
  conversion can be disabled using \code{map.method=NA}, as
  it is by default when \code{x} is a \code{matrix}, whose
  rows are assumed to match the rows in the target matrix.}

  \item{method}{method to use to deconvolve the target data
  and estimate cell proportions. The method must be a
  deconvolution algorithm that is able to run using
  signatures as only auxiliary input.

  The default method is \sQuote{lsfit}, which implements
  the algorithm proposed by \cite{Abbas2009} that is based
  on standard regression. An alternative method is the
  quadratic programming approach from \cite{Gong2011},
  which solves a nonnegative least-square problem with
  sum-up-to one constraints on the proportions.}

  \item{subset}{optional subset of features to use in the
  estimation.}

  \item{map.method}{method used to convert the basis
  signature's identifiers to match the target data's own
  type of identifiers. See \code{\link{mapIDs}}. Identifier
  conversion can be disabled using \code{map.method=NA}.}

  \item{...}{extra arguments passed to \code{\link{ged}}}

  \item{log}{logical that indicates if the computation
  should take place in log or linear scale. If \code{TRUE},
  all non-log-scaled data (signatures and/or target) are
  log-transformed using with base \code{lbase}. If
  \code{FALSE}, all log-scaled data (signatures and/or
  target) are exp-transformed using with base \code{lbase}.
  If a number, then the function acts as if \code{log=TRUE}
  using the value of \code{log} as \code{lbase}. If
  \code{NULL}, then log-transform is applied only if either
  the signatures or the target data is in log scale,
  otherwise non-log-scaled data is exp-transformed into
  linear values, via \code{\link{expb}(A, lbase)}. If
  \code{log=NA} no transformation is performed at all.}

  \item{lbase}{numeric base use for the
  logarithmic/exponential transformations that are applied
  to the signature or data matrix.}

  \item{normalize}{character string that specifies the
  normalisation method to apply jointly to the combined
  data (signatures,data). The normalisation is performed
  after transforming the data and/or signatures if
  necessary.}

  \item{verbose}{logical that toggles verbosity. A number
  (integer) can be passed to specify the verbosity level
  (the higher the more messages are output). Passing
  \code{verbose=Inf} toggles debug mode (all messages).
  Note that because it appears after \code{...} it must be
  fully named.}
}
\description{
  \code{gedProportions} implements a pre-processing
  pipeline for applying deconvolution methodologies that
  use a known set of cell type-specific signatures in order
  to estimate cell proportions in heterogeneous samples
  (e.g., \cite{Abbas et al. (2009)} or \cite{Gong et al.
  (2011)}).
}
\details{
  The actual estimation is performed via the
  \code{\link{ged}} interface, using a suitable
  deconvolution method.

  Before calling \code{\link{ged}}, the following
  pre-processing pipeline is applied to the data and/or the
  signature matrix: \describe{ \item{map}{ the gene
  identifiers of the signature matrix into identifiers in
  the target global expression matrix;} \item{subset}{
  signatures and data matrices to a common set of
  features;} \item{transform}{ signatures and data matrices
  to a common scale: linear or log; Log-scale is
  automatically detected using the same heuristic as
  \code{GEO2R}.} \item{normalise}{ jointly the signatures
  and data matrices using quantile normalisation.} }

  All steps are optional and can be disabled if needed (see
  argument details).
}
\references{
  Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z and Clark
  HF (2009). "Deconvolution of blood microarray data
  identifies cellular activation patterns in systemic lupus
  erythematosus." _PloS one_, *4*(7), pp. e6098. ISSN
  1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0006098>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/19568420>.

  Gong T, Hartmann N, Kohane IS, Brinkmann V, Staedtler F,
  Letzkus M, Bongiovanni S and Szustakowski JD (2011).
  "Optimal deconvolution of transcriptional profiling data
  using quadratic programming with application to complex
  clinical blood samples." _PloS one_, *6*(11), pp. e27156.
  ISSN 1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0027156>, <URL:
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3217948\&tool=pmcentrez\&rendertype=abstract>.
}
\seealso{
  \code{\link{ged}}, \code{\link{gedBlood}}
}

