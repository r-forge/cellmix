\name{match.nmf}
\alias{match.cell}
\alias{match.nmf}
\title{Assigning Estimated Signature to Real Cell-Types}
\usage{
  match.nmf(object, reference, n = NA)

  match.cell(object, reference, n = NA)
}
\arguments{
  \item{object}{an \code{NMF} object}

  \item{reference}{a list of markers, usually a
  \code{MarkerList} object.}

  \item{n}{the number of markers to use per cell-type to
  actually assign the components}
}
\value{
  \code{object} with the basis components reordered -- and
  renamed -- according to the computed assigned cell-types.
}
\description{
  Assigns and reorder basis components to a cell type based
  on the proportion of cell-specific markers that
  contribute the most to each component.

  \code{match.cell} is an alias for \code{match.nmf}
}

