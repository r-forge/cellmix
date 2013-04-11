\name{gmarkers}
\alias{gMarkerList}
\alias{gmarkers}
\alias{rMarkerList}
\alias{rmarkers}
\title{Generating Marker Lists}
\usage{
  gmarkers(n, k = -3, names = FALSE, values = FALSE,
    shuffle = FALSE)

  gMarkerList(...)

  rmarkers(...)

  rMarkerList(...)
}
\arguments{
  \item{n}{number of marker sets to generate, e.g., the
  number of cell types in the data.}

  \item{k}{numeric (integer) value that specifies the
  number of markers per set. If a negative value is
  specified, then it is interpreted as a growing factor for
  successive sets of markers, e.g., the default \code{k=-3}
  generates a first set of 3 markers (=3 x 1), a second set
  with 6 (= 2 x 3), a third with 9 (= 3 x 3), and so long.}

  \item{names}{specification for names. If a logical it
  indicates if names should be generated or not. If a
  single string, it indicates the prefix to use for the
  names, the default being \dQuote{Marker.}. If a character
  vector, then it specifies the names of each marker, and
  must then contain at least as many unique elements as
  markers in the whole generated list.}

  \item{values}{specification for values. If a logical, it
  indicates if a value should be generated and associated
  to each marker. If a numeric vector, then it specifies
  the values to associate with each marker. The names of
  this vector are currently not used. If not shuffling (see
  argument \code{shuffle}) and \code{values} is a numeric
  vector, then it must contain at least as many unique
  elements as markers in the whole generated list.}

  \item{shuffle}{logical that indicate if the markers
  should be shuffled randomly. This is used only if fixed
  names or values have been specified.}

  \item{...}{arguments eventually passed to
  \code{gmarkers}.}
}
\description{
  The functions documented here serve to generate marker
  lists and \code{\linkS4class{MarkerList}} objects. They
  are useful to develop and test deconvolution methods that
  use marker gene lists.

  \code{gmarkers} is the work horse generator function. It
  can generate "fixed" or random plain list of marker,
  optionnally with names and/or values.

  This is a shortcut for \code{MarkerList(gmarkers(...))}

  It is a shorcut for \code{gmarkers(..., shuffle=TRUE)}.

  This is a shortcut for \code{MarkerList(rmarkers(...))}
}

