\name{subsetML}
\alias{subsetML}
\title{Subsetting Data with MarkerList Objects}
\usage{
  subsetML(x, i, j, ..., drop = FALSE)
}
\arguments{
  \item{drop}{logical passed to the subset method
  \code{'['} of \code{x}.}

  \item{x}{ object from which to extract element(s) or in
  which to replace element(s).  }

  \item{i}{ indices specifying elements to extract or
  replace.  Indices are \code{numeric} or \code{character}
  vectors or empty (missing) or \code{NULL}.  Numeric
  values are coerced to integer as by
  \code{\link{as.integer}} (and hence truncated towards
  zero).  Character vectors will be matched to the
  \code{\link{names}} of the object (or for
  matrices/arrays, the \code{\link{dimnames}}): see
  \sQuote{Character indices} below for further details.

  For \code{[}-indexing only: \code{i}, \code{j},
  \code{\dots} can be logical vectors, indicating
  elements/slices to select.  Such vectors are recycled if
  necessary to match the corresponding extent.  \code{i},
  \code{j}, \code{\dots} can also be negative integers,
  indicating elements/slices to leave out of the selection.

  When indexing arrays by \code{[} a single argument
  \code{i} can be a matrix with as many columns as there
  are dimensions of \code{x}; the result is then a vector
  with elements corresponding to the sets of indices in
  each row of \code{i}.

  An index value of \code{NULL} is treated as if it were
  \code{integer(0)}.  }

  \item{j}{ indices specifying elements to extract or
  replace.  Indices are \code{numeric} or \code{character}
  vectors or empty (missing) or \code{NULL}.  Numeric
  values are coerced to integer as by
  \code{\link{as.integer}} (and hence truncated towards
  zero).  Character vectors will be matched to the
  \code{\link{names}} of the object (or for
  matrices/arrays, the \code{\link{dimnames}}): see
  \sQuote{Character indices} below for further details.

  For \code{[}-indexing only: \code{i}, \code{j},
  \code{\dots} can be logical vectors, indicating
  elements/slices to select.  Such vectors are recycled if
  necessary to match the corresponding extent.  \code{i},
  \code{j}, \code{\dots} can also be negative integers,
  indicating elements/slices to leave out of the selection.

  When indexing arrays by \code{[} a single argument
  \code{i} can be a matrix with as many columns as there
  are dimensions of \code{x}; the result is then a vector
  with elements corresponding to the sets of indices in
  each row of \code{i}.

  An index value of \code{NULL} is treated as if it were
  \code{integer(0)}.  }

  \item{...}{ indices specifying elements to extract or
  replace.  Indices are \code{numeric} or \code{character}
  vectors or empty (missing) or \code{NULL}.  Numeric
  values are coerced to integer as by
  \code{\link{as.integer}} (and hence truncated towards
  zero).  Character vectors will be matched to the
  \code{\link{names}} of the object (or for
  matrices/arrays, the \code{\link{dimnames}}): see
  \sQuote{Character indices} below for further details.

  For \code{[}-indexing only: \code{i}, \code{j},
  \code{\dots} can be logical vectors, indicating
  elements/slices to select.  Such vectors are recycled if
  necessary to match the corresponding extent.  \code{i},
  \code{j}, \code{\dots} can also be negative integers,
  indicating elements/slices to leave out of the selection.

  When indexing arrays by \code{[} a single argument
  \code{i} can be a matrix with as many columns as there
  are dimensions of \code{x}; the result is then a vector
  with elements corresponding to the sets of indices in
  each row of \code{i}.

  An index value of \code{NULL} is treated as if it were
  \code{integer(0)}.  }
}
\description{
  The function \code{subsetML} enables direct subsetting of
  data objects with \code{\link{MarkerList}} objects, as
  long as they have a suitable \code{\link{featureNames}},
  \code{\link{rownames}} or \code{\link{names}} method, as
  well as a subset method \code{[}. It is used for defining
  \code{'['} methods for \code{\linkS4class{NMF}} and
  \code{\linkS4class{ExpressionMix}}.
}
\examples{
data(sample.ExpressionSet)
e <- sample.ExpressionSet
annotation(e)
}

