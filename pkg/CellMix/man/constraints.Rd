\name{maxConstraint}
\alias{maxConstraint}
\title{Building Constraint Matrices}
\usage{
  maxConstraint(n, i_max = 1, alpha = 1)
}
\arguments{
  \item{n}{number of components.}

  \item{i_max}{index of the component that must be greater
  that all others.}

  \item{alpha}{minimum ratio between the \code{k}-th and
  all the other components.}
}
\description{
  These are utility functions to work with constraint
  matrices, which define linear equality and inequality
  constraints in optimization problems.

  \code{maxConstraint} creates a constraint matrix that
  specifies that a given component must be greater than all
  the others. It is used, in particular, in the
  deconvolution method from \cite{Gong et al. (2011)} to
  enforce which is\dQuote{qprog}
}
\details{
  \code{maxConstraint} creates a constraint matrix that
  specifies that a component \eqn{k} must be greater than
  all the others:

  \deqn{x_k \geq \alpha x_i }{x_k >= alpha * x_i},

  for all \eqn{1\leq i\leq n}{1<= i <= n}, where
  \eqn{\alpha}{alpha} is a fixed parameter is the minimum
  ratio between the \eqn{k}-th and all the others
  components, i.e. it controls the degree of differences
  between these components.
}
\examples{
maxConstraint(5, 3)
maxConstraint(5, 3, alpha=3)
}
\references{
  Gong T, Hartmann N, Kohane IS, Brinkmann V, Staedtler F,
  Letzkus M, Bongiovanni S and Szustakowski JD (2011).
  "Optimal deconvolution of transcriptional profiling data
  using quadratic programming with application to complex
  clinical blood samples." _PloS one_, *6*(11), pp. e27156.
  ISSN 1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0027156>, <URL:
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3217948\&tool=pmcentrez\&rendertype=abstract>.
}
\keyword{internal}

