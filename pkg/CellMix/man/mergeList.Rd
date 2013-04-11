\name{mergeList}
\alias{mergeList}
\title{Merge two lists}
\usage{
  mergeList(x, y, overwrite = FALSE)
}
\arguments{
  \item{x}{list to which elements from \code{y} are added.}

  \item{y}{list whose elements are added to \code{x}}

  \item{overwrite}{a logical indicating whether the
  elements of \code{y} already present in \code{x} should
  overwrite the latter or not (default).}
}
\description{
  Adds elements from a list to another list, optionally
  overwriting them if already present.
}

