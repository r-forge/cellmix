\name{DSAproportions}
\alias{DSAproportions}
\alias{.DSAproportions}
\title{Digital Sorting Algorithm: Proportion Estimation Method}
\usage{
  DSAproportions(x, markers, verbose = FALSE)

  .DSAproportions(x, markers, verbose = FALSE)
}
\arguments{
  \item{x}{target mixed expression matrix}

  \item{markers}{marker specification}

  \item{verbose}{verbose level}
}
\value{
  a matrix of proportions with the same number of columns
  as \code{x} and as many rows as elements in
  \code{markers} (i.e. cell types).
}
\description{
  \code{DSAproportions} implements the method proposed by
  \cite{Zhong et al. (2013)} as part of the
  \code{\link[=DSA-ged]{DSA}} algorithm. This method
  estimates cell proportions from mixed sample expression
  data, given a set of markers, i.e. features that are
  known to be exclusively expressed by a single cell type
  in the mixture.

  \code{.DSAproportions} is the internal -- non-exported --
  function that actually performs the estimation. It
  expects a list of valid indexes for each cell type whose
  proportions are to be estimated.
}
\examples{
# random data
x <- rmix(3, 100, 20)
m <- getMarkers(x)

# estimate proportions
p <- DSAproportions(x, m)
# plot against true proportions
profplot(x, p)
}
\references{
  Zhong Y, Wan Y, Pang K, Chow LM and Liu Z (2013).
  "Digital sorting of complex tissues for cell
  type-specific gene expression profiles." _BMC
  Bioinformatics_, *14*(1), pp. 89. ISSN 1471-2105, <URL:
  http://dx.doi.org/10.1186/1471-2105-14-89>, <URL:
  http://www.biomedcentral.com/1471-2105/14/89>.
}

