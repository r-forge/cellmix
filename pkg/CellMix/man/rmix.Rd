\name{rmix}
\alias{rmix}
\title{Generating Random Global Mixed Gene Expression Data}
\usage{
  rmix(x, n = 100, p = 20, markers = ceiling(nrow(x)/20),
    min = 0, max = 20, mfold = 2, alpha = 1,
    snoise = list(mean = 0, sd = 0.05),
    gnoise = list(shape = 5, scale = 1), ...)
}
\arguments{
  \item{x}{number of true underlying cell types or a matrix
  containing the signatures themselves, i.e. cell-specific
  expression values for each feature. For convenience, it
  may also specify the markers to enforce on the
  signatures, as a vector or list of length > 1, in which
  case argument \code{markers} must be missing.}

  \item{n}{number of features, i.e. genes. The argument is
  required if \code{x} specifies the number of signatures.
  If \code{x} is provided as a matrix, then \code{n} is
  used to subset it (\code{x[n, ]}) before simulating the
  global expression data.}

  \item{p}{number of samples}

  \item{markers}{specification of the number of markers to
  enforce on each cell type signature. This should be a
  value supported by \code{\link{enforceMarkers}}. Markers
  enforcement may be disabled with \code{markers=NA}.}

  \item{min}{minimum cell-specific expression value before
  adding noise and marker differential.}

  \item{max}{maximum cell-specific expression value before
  adding noise and marker differential.}

  \item{mfold}{fold change expected on cell-specific
  expression for marker genes}

  \item{alpha}{parameter for the dirichlet distribution
  from which are drawn the mixture proportions, using
  \code{\link[gtools]{rdirichlet}}.}

  \item{snoise}{parameters for the normal noise added to
  each true underlying signatures as \eqn{x + N(\mu,
  \sigma)}.}

  \item{gnoise}{parameters for the normal noise with
  inverse gamma variance added to each feature global
  expression profile as \eqn{e_{ij} + N(0, 1/\gamma_i)}.}

  \item{...}{extra arguments currently not used.}
}
\value{
  an \code{\linkS4class{ExpressionMix}} object, that
  contains the true underlying signatures and proportions
  stored as an NMF model.
}
\description{
  The function \code{rmix} generates an
  \code{\linkS4class{ExpressionMix}} object, composed of a
  given number of underlying cell types. The amount of
  noise added to both the cell-specific signatures and the
  global expression values is customisable.
}
\examples{
# 3 cell types, 100 features, 20 samples
rmix(3, 100, 20)

# from known signature matrix
s <- rmatrix(100, 5)
x <- rmix(s, p=20)
dim(x)
if( !isCRAN_timing() ){
aheatmap(x)
}

# markers are enforced on each true signature
x <- rmix(4, 50, 20, markers=6)
if( !isCRAN_timing() ){
basismap(x, Rowv=NA)
}
# or also
x <- rmix(1:4, 50, 20)
if( !isCRAN_timing() ){
basismap(x, Rowv=NA)
}
}

