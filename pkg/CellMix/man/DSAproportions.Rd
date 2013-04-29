\name{DSAproportions}
\alias{DSAproportions}
\alias{.DSAproportions}
\title{Digital Sorting Algorithm: Proportion Estimation Method}
\usage{
  DSAproportions(x, markers, log = NULL,
    match.names = TRUE, verbose = FALSE)

  .DSAproportions(x, markers)
}
\arguments{
  \item{x}{target mixed expression matrix}

  \item{markers}{marker specification}

  \item{log}{indicates if the data are in log-scale or
  should be converted to linear-scale. This is relevant
  because the DSA algorithm assumes that the input mixed
  data are in linear scale (i.e. not log-trasnformed).

  If \code{NULL}, then data's scale is detected by
  \code{link{is_logscale}} and conversion to linear-scale
  is performed if necessary. If \code{TRUE} the data is
  exponentialised (using log base 2). If \code{FALSE} the
  data is left unchanged (the detected log scale is
  displayed in verbose mode). If a number, then it is used
  as the log base to exponentialise the data.}

  \item{match.names}{logical that indicates if marker names
  should be first matched against the data, to convert them
  into integer indexes. If code{FALSE}, then \code{markers}
  should be a list of integer vectors that corresponds to
  the indexes of markers for each cell types.}

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

