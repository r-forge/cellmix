\docType{methods}
\name{MarkerList}
\alias{MarkerList}
\alias{MarkerList,ANY-method}
\alias{MarkerList,character-method}
\alias{MarkerList,ExpressionSet-method}
\alias{MarkerList,factor-method}
\alias{MarkerList,integer-method}
\alias{MarkerList,list-method}
\alias{MarkerList,MarkerList-method}
\alias{MarkerList,matrix-method}
\alias{MarkerList-methods}
\alias{MarkerList,missing-method}
\alias{MarkerList,vector-method}
\title{Factory Method for Marker Lists}
\usage{
  MarkerList(object, ...)

  \S4method{MarkerList}{missing}(object, ..., file = NULL)

  \S4method{MarkerList}{vector}(object, ..., names = NULL)

  \S4method{MarkerList}{character}(object, ...,
    names = NULL)

  \S4method{MarkerList}{list}(object, ..., unlist = FALSE,
    quiet = FALSE)

  \S4method{MarkerList}{matrix}(object, ..., values = TRUE,
    add.names = TRUE)
}
\arguments{
  \item{object}{input object from which marker data are
  extracted.}

  \item{...}{extra arguments to allow extension. See each
  method's description for more details.}

  \item{file}{filename of a file that contains marker data.
  Conneciton objects are also supported.}

  \item{names}{factor or character vector that defines the
  types associated to each marker in \code{object}.}

  \item{unlist}{logical that indicates if the result should
  be a vector representation of the marker list. In this
  case auxiliary numeric values are dropped from the list.}

  \item{quiet}{logical that indicates if warnings should be
  silenced if any.}

  \item{values}{a logical that indicates if the row
  maximums should be embedded into the marker list.}

  \item{add.names}{logical that indicates if names should
  be added if necessary (\code{TRUE}). This is used in the
  case the input matrix has no rownames and the row maximum
  values are requested to be attached to the marker list
  (\code{values=TRUE}), which throws an error if
  \code{add.names=FALSE}.}
}
\description{
  \code{MarkerList} is an S4 generic function that provides
  a convenient interface to create
  \code{\linkS4class{MarkerList}} objects from a variety of
  input formats, such as plain named lists, factors,
  matrices or \code{\linkS4class{ExpressionSet}} objects,
  text files, etc..
}
\details{
  Except for its method \code{MarkerList,matrix}, this
  function does not try to infer markers from the input
  data, but only reorganise groups of markers that are
  already clearly defined in the input data.
}
\section{Methods}{
  \describe{

  \item{MarkerList}{\code{signature(object =
  "MarkerList")}: This method acts as a Constructor-Copy
  that returns the input MarkerList object unchanged if no
  other argument are passed, or creates a new MarkerList
  object based on \code{object}, using the arguments in
  \code{...} to initialise the slots. }

  \item{MarkerList}{\code{signature(object = "missing")}:
  Creates a \code{MarkerList} object as per
  \code{new('MarkerList', ...)} or from a file (see section
  \emph{Loading from files}). }

  \item{MarkerList}{\code{signature(object = "ANY")}:
  Default method that tries to retrieve marker data from
  the input object using the function
  \code{\link{getMarkers}}. }

  \item{MarkerList}{\code{signature(object = "vector")}:
  Creates a \code{MarkerList} object treating the input
  vector as a factor, except if the input is a numeric
  vector of scores. See
  \code{\link{MarkerList,factor-method}}.

  If \code{object} is numeric, one also needs to supply a
  factor that defines the cell types. Markers ids are taken
  from the names \code{object}, or, if missing, from the
  names of \code{types}. }

  \item{MarkerList}{\code{signature(object = "integer")}:
  Creates a MarkerList object using the -- normally --
  duplicated names of the input vector as types, and the
  values as marker indexes. See
  \code{\link{MarkerList,factor-method}}. }

  \item{MarkerList}{\code{signature(object = "character")}:
  Either loads a marker list from the internal registry
  (see \code{\link{cellMarkers}}), or use the names of
  \code{object} for defining the sets and the values
  (strings) as marker identifiers.

  If no names are present, then it uses the character
  vector as a factor, whose levels are the sets, and
  creates an index marker list with
  \code{MarkerList,factor-method}. }

  \item{MarkerList}{\code{signature(object = "factor")}:
  Create a MarkerList object by splitting the input factor
  into a list of its different levels.

  If \code{object} has names they are used as marker
  identifiers, otherwise an index marker list is created. }

  \item{MarkerList}{\code{signature(object = "list")}:
  Convert a standard list into a MarkerList object.

  The list should contain no duplicated marker identifiers
  (either as strings or integer indexes). }

  \item{MarkerList}{\code{signature(object = "matrix")}:
  Creates a MarkerList object from a numeric matrix, that
  is assumed to have one column per marker set, e.g., a
  matrix cell type-specific signature.

  This method uses a very simple heuristic, which associate
  each row in a matrix to the column with the maximum entry
  (see \code{\link{max.col}}).

  For more complex and finer ways of extracting markers
  from expression data see \code{\link{extractMarkers}}.

  The rownames are used as marker identifiers if present.
  If not otherwise specified, each marker is associated
  with its correspondong row maximum, which is also stored
  in the result object. }

  \item{MarkerList}{\code{signature(object =
  "ExpressionSet")}: Create a MarkerList from the
  expression matrix of an
  \code{\linkS4class{ExpressionSet}} object.

  This is a shortcut for \code{MarkerList(exprs(object),
  ...)}. }

  }
}
\section{Loading from files}{
  Loading marker lists directly from files can be performed
  by providing filename to argument \code{file}. Currently
  this functionality is only implemented for plain text
  files that contain one marker feature description per
  line: marker names (e.g., ENTREZID), marker set (e.g.,
  cell type name) and optionally some numeric value (e.g.
  specificity score).

  In this case all arguments in \code{...} are passed to
  \code{\link{read.table}}, and should be set so that the
  loaded data.frame conforms with one of the following
  formats: \describe{ \item{at least one column and row
  names:}{ the row names are used as marker names, column 1
  as marker set. Note that tables with automatic rownames
  generated by \code{read.table} (i.e. 1,2,3, etc..) are
  not considered as proper rownames.} \item{two columns and
  no row names:}{column 1 is used for marker names, column
  2 for the marker set names.} \item{single column and no
  (or automatic) rownames:}{column is used for marker set
  names, the markers are taken to be integers.} }
  Optionally, in all cases, the next numeric column after
  the one used for marker set names is used for associated
  numeric values.
}
\examples{
## create a MarkerList from another MarkerList
m <- rMarkerList(3, 4, names=TRUE)
m
MarkerList(m, setName=mkScalar("Another list of markers"))
txtdesc <- function(x) textConnection(paste(x, collapse="\\n"))

# two columns
ml <- MarkerList(file=txtdesc(c("m1 a", "m2 b", "m3 a")))
summary(ml)

# single column with row.names
ml2 <- MarkerList(file=txtdesc(c("m1 a", "m2 b", "m3 a")), row.names=1L)
identical(ml, ml2)

# three columns with values
mlv <- MarkerList(file=txtdesc(c("m1 a 1.4", "m2 b 1.5", "m3 a 2.6")))
summary(mlv)
print(mlv)

# skip an extra column
mlv <- MarkerList(file=txtdesc(c("m1 a 1.4", "m2 b 1.5", "m3 a 2.6")))
summary(mlv)
print(mlv)

# single column: index MarkerList
ml <- MarkerList(file=txtdesc(c("a", "b", "a")))
summary(ml)
## create a MarkerList from an integer vector
# names are required
try( MarkerList( 1:15 ) )
# repeated names define marker sets
MarkerList( setNames(1:12, rep(letters[1:3], 4)) )
## create a MarkerList from a character vector
# no names: repeated values define the sets
v <- sample(letters[1:3], 15, replace=TRUE)
MarkerList(v)
# with names: repeated names define the sets
m <- str_c('M', 1:15)
MarkerList( setNames(m, v) )
## create a MarkerList from a factor
f <- factor(sample(letters[1:3], 15, replace=TRUE))
MarkerList(f)
MarkerList( addNames(f, 'M') )
## create an index MarkerList from a standard list
m <- list(1:3, 5:10, c(4L,12L))
MarkerList(m)

# with mixed integer/numeric, a warning is thrown
m <- list(1:3, 5:10, c(4,12))
MarkerList(m)
## create a MarkerList from a matrix of signatures
s <- rmatrix(10,4)
m <- MarkerList( s )
# the object has values, that can be dropped is necessary
m
dropvalues(m)
# or directly
identical(dropvalues(m), MarkerList(s, values=FALSE))
}
\keyword{methods}

