\name{hasAnnotation}
\alias{getAnnotation}
\alias{hasAnnotation}
\title{Extracting Annotation from Objects}
\usage{
  hasAnnotation(object, ...)

  getAnnotation(object, ..., null = TRUE)
}
\arguments{
  \item{object}{an object}

  \item{...}{extra parameters (currently not used)}

  \item{null}{logical that indicates if an empty character
  string should be return as \code{NULL}.}
}
\description{
  \code{hasAnnotation} tells if an object has some -- non
  empty -- attached annotation.

  \code{getAnnotation} try extracting embedded annotations.
  By default, it returns \code{NULL} if the object contains
  no annotation.
}

