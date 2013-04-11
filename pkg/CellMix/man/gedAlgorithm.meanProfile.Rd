\name{gedAlgorithm.meanProfile}
\alias{gedAlgorithm.meanProfile}
\alias{meanProfile-ged}
\title{Partial Gene Expression Deconvolution: Marker Mean Expression Profile}
\usage{
  gedAlgorithm.meanProfile(y, x, scale = FALSE, ...)
}
\arguments{
  \item{y}{target expression matrix}

  \item{x}{specification of which feature to use as
  markers, as a list, typically a \code{\link{MarkerList}}
  object.}

  \item{scale}{logical that indicates if the proportion
  estimates should be scale to sum-up to one before
  returning them.}

  \item{...}{extra arguments not used.}
}
\description{
  The algorithm \sQuote{meanProfile} uses a set of known
  marker genes for each cell type to compute a mean
  expression profile within each cell type separately.

  This is in essence the preliminary step proposed by
  \cite{Kuhn et al. (2011)} in order to compute a proxy for
  the actual cell proportions.

  The average profiles are expected to be correlate well
  with the actual proportions, provided that the individual
  markers gene expression profiles are not too noisy and
  that the markers are indeed markers.
}
\details{
  \strong{Important:} this method does not compute
  cell-specific differential expression as described in
  \cite{Kuhn et al. (2011)}, but only the cell proportion
  proxy. Hence, the result only contains estimated
  proportions (accessible via \code{\link{coef}}) and an
  empty basis signature matrix.
}
\examples{
# random data with markers
x <- rmix(3)
m <- getMarkers(x)

# compute proxy proportions
res <- ged(x, m, method='meanProfile')
# no cell-specific signatures
dim(res)

#NB: estimates are not scaled to sum up to one
profplot(x, res)
profplot(x, res, scale=TRUE)
}
\references{
  Kuhn A, Thu D, Waldvogel HJ, Faull RLM and Luthi-Carter R
  (2011). "Population-specific expression analysis (PSEA)
  reveals molecular changes in diseased brain." _Nature
  methods_, *8*(11), pp. 945-7. ISSN 1548-7105, <URL:
  http://dx.doi.org/10.1038/nmeth.1710>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/21983921>.
}

