\docType{class}
\name{MarkerList-class}
\alias{attachMarkers}
\alias{getMarkers}
\alias{has.markers}
\alias{isMarkerList}
\alias{MarkerList-class}
\alias{summary,MarkerList-method}
\title{Class for Marker Gene Lists}
\usage{
  isMarkerList(x)

  \S4method{summary}{MarkerList}(object)

  attachMarkers(x, object)

  getMarkers(x, err.notfind = NULL, unlist = FALSE)

  has.markers(x)
}
\arguments{
  \item{x}{generally a \code{MarkerList} object or any
  object depending on the method. See their respective
  description for details.}

  \item{object}{a list of markers, usually a
  \code{MarkerList} object, or any object with an attached
  list of markers (cf. \code{attachMarkers} and
  \code{getMarkers}), depending on the method. See their
  respective description for details.}

  \item{err.notfind}{a character vector that will prefix an
  error message if no \code{x} has no marker list attached.
  The error is only thrown if \code{err.notfind} is not
  \code{NULL}.}

  \item{unlist}{a single logical specifying if the result
  should be unlisted Default is \code{FALSE}.}
}
\value{
  \item{summary}{ returns a numeric vector containing the
  number of markers for each cell-type}

  \item{\code{attachMarkers}}{returns \code{object} with an
  attribute named \code{'markers'} that contains the result
  of \code{MarkerList(x)}.}

  \item{\code{getMarkers}}{returns the attached marker
  list}

  \item{\code{has.markers}}{returns a logical}
}
\description{
  This is the main class used to store sets of marker
  genes, for use in gene expression deconvolution methods.
  It serves as a light version of class
  \code{\linkS4class{MarkerSetCollection}}, assuming that
  all marker identifiers are of the same type (e.g., entrez
  gene, unigene or probe ids), and are relative to the same
  platform, organism, etc.. (see details).

  \code{isMarkerList} tests if an object is a
  \code{\linkS4class{MarkerList}} object.

  \code{attachMarkers} attaches a marker list to an
  \code{NMF} object, by storing it in an attribute named
  \code{'markers'}. It actually stores the results from
  \code{MarkerList(x)}, that is a \code{MarkerList} object.

  \code{getMarkers} extracts the marker list from the
  attribute named \code{'markers'} that will generally
  contain a \code{MarkerList} object.

  \code{has.markers} tells if an object has an attached
  marker list, i.e. an attribute named \code{'markers'}.
}
\details{
  The objective of the \code{MarkerList} class is to
  simplify the structure and processing of gene lists,
  compared to handling \code{MarkerSetCollection} objects.
  Conversion methods between these two classes are
  provided.

  The class is essentially a named list, in which each
  element contains data about a set of marker
  features/genes, e.g., from a given cell-type. Marker
  features are assumed to be exclusive to each set, i.e.
  they appear in only one element of the list. In addition
  to their identifiers, markers can be associated with
  numeric values, e.g., corresponding to p-values or
  specificity scores, or integer, which are then
  interpreted as indexes relative to some expression data.

  It contains the same slots as
  \code{\linkS4class{GeneSet}}, which are used when
  converting \code{MarkerList} objects into
  \code{MarkerSetCollection} objects, to fill the slots of
  all gene set. For a description of each slot, please see
  the documentation for \code{\linkS4class{GeneSet}}.
}
\section{Methods}{
  \describe{

  \item{[}{\code{signature(x = "MarkerList", i = "ANY", j =
  "ANY")}: \code{MarkerList} objects have convenient subset
  methods, that allow to subset the list not only along the
  sets (first dimension) but also along the markers, which
  proves to be very useful to subset markers list based on
  another set of identifier, e.g., the genes present in
  some expression data, or in another marker list.

  The method \code{`[`} works performs basic strict subsets
  using integer, character and logical subsetting vectors.

  See \code{\link{[,MarkerList,ANY,ANY-method}} for more
  details. }

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

  \item{annotation}{\code{signature(object =
  "MarkerList")}: Extract name of the annotation package
  from a marker list. }

  \item{annotation<-}{\code{signature(object =
  "MarkerList", value = "character")}: Sets the name of the
  annotation package associated to a marker list. }

  \item{annotation<-}{\code{signature(object =
  "MarkerList", value = "NULL")}: Reset the annotation
  package associated to a marker list }

  \item{asCBC}{\code{signature(object = "MarkerList")}:
  This method combines markers of cell types that belong to
  the same CBC group. }

  \item{.atrack}{\code{signature(object = "MarkerList")}:
  This method allows to add annotation tracks in heatmaps
  produced by \code{\link{aheatmap}}, to highlight the
  position of markers, when plotting either the global
  expression values or cell type-specific signatures,
  whether measured or estimated. }

  \item{combine}{\code{signature(x = "MarkerList", y =
  "factor")}: Combine markers from multiple cell types of a
  MarkerList object, based on groups defined by a given
  factor. }

  \item{Compare}{\code{signature(e1 = "MarkerList", e2 =
  "numeric")}: \code{Compare} compares all embedded values
  with a given fixed value. This is useful to filter
  markers based on their associated scores. }

  \item{connectivity}{\code{signature(object =
  "MarkerList")}: Computes the connectivity matrix of each
  set (i.e. cell/tissue type) of a marker list, i.e. the
  matrix whose entries are the number of features in common
  between each pair o cell type. }

  \item{convertIDs}{\code{signature(object = "MarkerList",
  to = "GeneIdentifierType", from = "GeneIdentifierType")}:
  Convert IDs from a MarkerList object.

  In this case, argument \code{unlist} indicates if the
  result should be a simple list containing the mapping (a
  list) for each cell type or a
  \code{\linkS4class{MarkerList}} object (default). }

  \item{details}{\code{signature(object = "MarkerList")}:
  Shows more details than the regular \code{show} method. }

  \item{drop}{\code{signature(x = "MarkerList")}:
  \code{drop} drops empty sets of markers from a
  \code{MarkerList} object, as well as markers with
  \code{FALSE} values if the object is a logical marker
  list. It returns the reduced \code{MarkerList} object. }

  \item{featureNames}{\code{signature(object =
  "MarkerList")}: Returns the character vector of all
  marker names in the list. }

  \item{flatten}{\code{signature(object = "MarkerList")}:
  For \code{\link{MarkerList}} objects the conversion is
  similar to what \code{\link{unlist}} would do, but
  argument \code{use.names} is used slightly differently.

  Secondly, the names of the main elements are append to
  the marker names (like \code{unlist} does if
  \code{use.names=TRUE}) only for marker lists that contain
  numeric values (e.g. specificity scores). In the other
  cases, they are used either as values for character
  marker lists or as duplicated names for integer marker
  lists. }

  \item{geneIdType}{\code{signature(object =
  "MarkerList")}: Returns the type of gene identifier used
  in a marker list. }

  \item{geneIdType<-}{\code{signature(object =
  "MarkerList", value = "NULL")}: Sets the type of gene
  identifier used in a marker list to
  \code{\link{NullIdentifier}()}. }

  \item{geneIdType<-}{\code{signature(object =
  "MarkerList", value = "character")}: Sets the type of
  gene identifier used in a marker list to a given value. }

  \item{geneIdType<-}{\code{signature(object =
  "MarkerList", value = "GeneIdentifierType")}: Sets the
  type of gene identifier associated with a marker list to
  a given \code{\linkS4class{GeneIdentifierType}} object. }

  \item{geneValues}{\code{signature(object =
  "MarkerList")}: Returns the values embedded in a
  \code{MarkerList} object. }

  \item{incidence}{\code{signature(x = "MarkerList")}:
  Computes the incidence matrix of each feature in each set
  (i.e. cell/tissue type) of a marker list.

  The matrix is computed by converting the
  \code{MarkerList} object into a
  \code{\link{GeneSetCollection}} object, for which a
  method \code{\link{incidence}} exists in the
  \pkg{GSEABase} package. }

  \item{MarkerList}{\code{signature(object =
  "MarkerList")}: This method acts as a Constructor-Copy
  that returns the input MarkerList object unchanged if no
  other argument are passed, or creates a new MarkerList
  object based on \code{object}, using the arguments in
  \code{...} to initialise the slots. }

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "matrix")}: Workhorse method is for
  \code{markermap}. }

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "ExpressionSet")}: The method \code{markermap} for
  \code{ExpressionSet} objects calls the main
  \code{markermap} method on the expression matrix
  \code{exprs(object)}. }

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "NMF")}: The method \code{markermap} for
  \code{NMF} objects calls the main \code{markermap} method
  on the basis matrix \code{basis(object)}. }

  \item{markermap}{\code{signature(object = "MarkerList",
  data = "NMFfitX")}: The method \code{markermap} for
  \code{NMFfitX} objects calls the main \code{markermap}
  method on the best fit \code{fit(object)}. }

  \item{reorder}{\code{signature(x = "MarkerList")}: The
  method \code{reorder} for \code{MarkerList} objects
  allows to reorder each set of markers according to some
  auxiliary reference data, e.g, their gene expression in
  different cell types.

  See \code{\link{reorder,MarkerList-method}} for more
  details. }

  \item{sapply}{\code{signature(X = "MarkerList")}: Applies
  a given function to a \code{\linkS4class{MarkerList}}
  object. }

  \item{subset}{\code{signature(x = "MarkerList")}:
  \code{subset} subset a \code{MarkerList} object keeping
  only the markers that are present with a given reference
  set, which can be a character vector, the rownames of a
  matrix or an \code{NMF} object, or the feature names of
  an \code{ExpressionSet}. The markers are matched using
  the function \code{matchIndex}, and gene identifier
  conversion is performed if necessary. }

  \item{summary}{\code{signature(object = "MarkerList")}:
  Summary method for \code{MarkerList} objects }

  \item{summary}{\code{signature(object = "MarkerList")}:
  Summary method for \code{MarkerList} objects }

  }
}
\section{Heterogeneous marker lists}{
  Due to the structure of \code{MarkerList} objects, all
  gene identifiers must share the same set of
  characteristics. However, it is possible to associate a
  \code{MarkerList} object with multiple annotation
  packages, which will be correctly dealt with if all the
  list's identifiers are found in only one of them. For
  example the registered maker list \sQuote{Abbas} is
  associated with annotation packages \code{hgu133a.db} and
  \code{hgu133b.db}.
}
\author{
  Renaud Gaujoux
}
\keyword{methods}

