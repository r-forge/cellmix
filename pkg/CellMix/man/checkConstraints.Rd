\name{checkConstraints}
\alias{checkConstraints}
\title{Checking Linear Constraints}
\usage{
  checkConstraints(x, A, b, margin = 2L,
    op = c(">=", "<=", "=="))
}
\arguments{
  \item{x}{matrix or vector to check.}

  \item{A}{left-hand matrix}

  \item{b}{right-hand vector}

  \item{margin}{which dimension should be checked (1L: rows
  or 2L: columns).}

  \item{op}{comparison operator to use (\code{'>=', '<=',
  '='}), indicating the type of constraints: inequality
  (greater than or lower than) or equality.}
}
\value{
  a logical vector of length \code{nrow(A)}.
}
\description{
  Checks if a given matrix or vector satisfies a set of
  linear constraints. The constraints are specified by a
  left-hand matrix \eqn{A} and a right-hand vector \eqn{b}
  such that a matrix or vector \eqn{x} satisfies the
  constraints if an only if \deqn{A \cdot x \geq b}{A x >=
  b}.
}

