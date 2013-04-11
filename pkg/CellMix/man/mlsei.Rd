\name{mlsei}
\alias{mlsei}
\title{Multivariate Least Squares with Equalities and Inequalities}
\usage{
  mlsei(A, B, ..., fulloutput = FALSE, verbose = TRUE)
}
\arguments{
  \item{A}{numeric matrix containing the coefficients of
  the quadratic function to be minimised, \eqn{||Ax-B||^2};
  if the columns of \code{A} have a names attribute, they
  will be used to label the output.}

  \item{B}{numeric matrix containing the right-hand side of
  the quadratic function to be minimised.}

  \item{fulloutput}{if \code{TRUE}, also returns the
  covariance matrix of the solution and the rank of the
  equality constraints -- only used if \code{type = 1} (see
  \code{\link{lsei}}).}

  \item{verbose}{logical that toggles verbose messages}

  \item{...}{extra arguments all passed to
  \code{\link{lsei}}.}
}
\value{
  If \code{fulloutput=TRUE}, a list whose elements are the
  results returned by \code{lsei} on each column of
  \code{B}, otherwise a list with the following elements:
  \item{X}{a matrix of dimension \code{ncol(B) x ncol(A)}
  that contain the estimated solution of the least-square
  problem.} \item{IsError}{a logical that indicates if an
  error occured} \item{type}{the string "lsei", such that
  how the solution was obtained can be traced}
  \item{solutionNorm}{a vector, whose elements are the sum
  of absolute values of residuals of equalities and
  violated inequalities, as computed for each column of
  \code{B}.} \item{solutionNorm}{a vector, whose elements
  are the values of the minimised quadratic function at the
  solution, i.e. the value of ||Ax-b||^2, as computed for
  each column of \code{B}.}
}
\description{
  This function provides a multivariate version of function
  \code{\link{lsei}} from the \code{\link{limSolve}}
  package. It solves the following optimisation problem for
  each column b of B
}
\details{
  \deqn{\min_x ||Ax - b||^2}{min ||Ax - b||^2}

  subject to the equality constraints: E x = f

  and/or inequality constraints: G x >= h

  It applies \code{lsei} on each column of a right-hand
  side matrix, and automatically transpose the input and
  result matrices back and forth if needed.
}
\examples{
# random left-hand quadratic coefficients
A <- rmatrix(20, 3)
r <- ncol(A)
# random right-hand target matrix
B <- rmatrix(20, 10)
# constrained least-squares: nonnegative and sum up to one
res <- mlsei(A,B, E=matrix(1, 1, r), F=1, G=diag(1, r), H=rep(0,r))
}

