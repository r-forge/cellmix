\docType{methods}
\name{[,MarkerList,ANY,ANY-method}
\alias{[,MarkerList,ANY,ANY-method}
\alias{[,MarkerList,list,ANY-method}
\alias{[,MarkerList,missing,list-method}
\alias{subset,MarkerList-method}
\title{Subsetting Marker Lists}
\usage{
  \S4method{[}{MarkerList,ANY,ANY}(x, i, j, ...,
    drop = FALSE)

  \S4method{[}{MarkerList,list,ANY}(x, i, j, ...,
    match.names = TRUE, drop = FALSE)

  \S4method{[}{MarkerList,missing,list}(x, i, j, ...,
    drop = FALSE)

  \S4method{subset}{MarkerList}(x, subset, select,
    filter = NULL, skip = 0, total = FALSE, invert = FALSE,
    fixed = TRUE, ignore.case = FALSE, verbose = FALSE,
    ...)
}
\arguments{
  \item{drop}{logical that indicates if one should drop
  empty cell types before returning the result -- as a
  \code{MarkerList}.}

  \item{match.names}{logical that indicates one should try
  to match the names of \code{x} and \code{i} before
  subsetting/re-ordering.}

  \item{subset}{the reference set in which the markers are
  looked-up. It can be a character vector, a matrix, or
  \code{NMF} or \code{ExpressionSet} objects.}

  \item{select}{a character, integer or logical vector
  specifying which set/cell-type to keep. If a single
  string and \code{fixed=FALSE}, it is used as a regular
  expression to match against the set names.}

  \item{filter}{a single integer specifying the maximum
  number of markers per cell-type}

  \item{skip}{a single integer specifying the number of top
  markers to skip in each cell type}

  \item{total}{logical that inidicates if argument
  \var{filter} specifies the total number of markers or the
  maximum number of markers per cell type.}

  \item{invert}{logical used when \code{select} that
  indicates that one selects the set (cell type) that do
  not satisfy the condition defined by argument
  \code{select}.}

  \item{fixed}{logical used only when \code{select} is a
  single string to specify that it should be interpreted as
  a regular expression to match against the set names.}

  \item{ignore.case}{logical used only when \code{select}
  is a single string to specify that the match should
  \strong{not} be case-sensitive.}

  \item{verbose}{verbosity level, usually \code{TRUE} or
  \code{FALSE}.}

  \item{...}{extra arguments passed to
  \code{\link{convertIDs}} if necessary.}

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
}
\value{
  an \code{MarkerList} object
}
\description{
  \code{MarkerList} objects have convenient subset methods,
  that allow to subset the list not only along the sets
  (first dimension) but also along the markers, which
  proves to be very useful to subset markers list based on
  another set of identifier, e.g., the genes present in
  some expression data, or in another marker list.

  The method \code{`[`} works performs basic strict subsets
  using integer, character and logical subsetting vectors.
}
\details{
  For \code{`[`}, argument \code{i} must be a
  \strong{valid} indexes, i.e. a character vector
  containing only strings that match exactly some names of
  \code{x}, or an integer with indexes in the range of the
  length of \code{x}, or a logical vector not longer than
  \code{x}.

  The second argument \code{j} is allowed not to match
  exactly, and may be used to subset the list against
  another set of identifiers, or limit each set to the top
  first markers.
}
\section{Methods}{
  \describe{

  \item{[}{\code{signature(x = "MarkerList", i = "list", j
  = "ANY")}: This method subsets/reorders each set of
  marker using a subset specification found in }

  \item{[}{\code{signature(x = "MarkerList", i = "missing",
  j = "list")}: This method is equivalent to \code{x[i, ,
  ..., match.names=FALSE]}, i.e. that each element of the
  marker list \code{x} are subset/reordered as its
  corresponding one in \code{j}, following the lists' order
  -- with no attempt to match their elements names before
  subsetting. }

  \item{subset}{\code{signature(x = "MarkerList")}:
  \code{subset} subset a \code{MarkerList} object keeping
  only the markers that are present with a given reference
  set, which can be a character vector, the rownames of a
  matrix or an \code{NMF} object, or the feature names of
  an \code{ExpressionSet}. The markers are matched using
  the function \code{matchIndex}, and gene identifier
  conversion is performed if necessary. }

  }
}
\examples{
# load IRIS markers and Abbas signature data
m <- MarkerList('IRIS')
m
# only keep the markers present in a dataset
data(Abbas)
m[, featureNames(Abbas)]
# NB: this could also be done with subset(m, Abbas)

# take top 10 (or smaller)
summary(m[,1:10])

# take markers with associated score >= 0.9
nmark(m[m >= 0.9])
}
\keyword{methods}

