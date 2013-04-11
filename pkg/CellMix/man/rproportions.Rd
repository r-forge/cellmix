\name{rproportions}
\alias{rproportions}
\title{Generating Random Cell Type Proportions}
\usage{
  rproportions(x, target, ..., alpha = 1)
}
\arguments{
  \item{...}{extra arguments passed to
  \code{\link[NMF]{rnmf}}.}

  \item{alpha}{shape parameter(s) for the dirichlet
  distribution. It must be a single numeric value, or a
  numeric vector of length the number of basis components
  in the NMF model -- which is generally specified in
  \var{object}. If a single numeric, it is used as a shape
  parameters for all components. If \code{alpha=NA}, then
  the coefficients are only scaled to sum-up to one.

  Note that because it appears after arguments \code{...},
  its name must be fully specified, or it will be passed to
  \code{rnmf}, while its default value used within
  \code{rproportions}.}

  \item{x}{an object that determines the rank, dimension
  and/or class of the generated NMF model, e.g. a numeric
  value or an object that inherits from class
  \code{\linkS4class{NMF}}. See the description of the
  specific methods for more details on the supported
  types.}

  \item{target}{optional specification of target
  dimensions. See section \emph{Methods} for how this
  parameter is used by the different methods.}
}
\description{
  \code{rproportions} generates a random NMF model, whose
  mixture coefficients for each sample (i.e. column of the
  coefficient matrix) fulfil the sum-up-to-one constraint
  of proportions.
}
\details{
  Internally, it generates a full random NMF model using
  the function \code{\link[NMF]{rnmf}}, and either scale
  the columns of its coefficient matrix so that they sum-up
  to one if \code{alpha=NA}, or re-draw them from a
  dirichlet distribution using the function
  \code{\link[gtools]{rdirichlet}} if \code{alpha} is a
  single number or a numeric vector.

  This function is suitable for seeding NMF computation,
  and is in fact registered as the seeding method
  \sQuote{rprop}, which can be used in calls to
  \code{\link[NMF]{nmf}}.
}
\examples{
# target data
x <- rmatrix(20,10)

# generate models
rproportions(3, x)
rproportions(3, 30, 15)

# scaling only
rproportions(3, 20, 10, alpha=NA)
# full shape specification
rproportions(3, x, alpha=c(1,2,3))

# error if alpha has wrong length
try( rproportions(3, 20, 10, alpha=c(1,2,3,4)) )
}

