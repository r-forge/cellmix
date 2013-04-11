\name{userData}
\alias{userData}
\title{User Data Directory}
\usage{
  userData(..., create = TRUE, package = "base")
}
\arguments{
  \item{...}{path parts passed to \code{\link{file.path}}
  to be appended to the main path.}

  \item{create}{logical that indicates if the directory
  should be created if it does not exists.}

  \item{package}{name of the package associated with the
  user data path. It is used to prefix the path, within the
  user R data directory.}
}
\description{
  \code{userData} returns the path to a local directory
  where package-related user data can be stored. Note that
  a directory is \strong{always} created if necessary (see
  details).
}
\details{
  If in interactive mode, the user is asked if the
  directory can be created in his home directory,
  otherwise, or if the user does not allow the creation in
  his home, the directory is created in the current R
  session's temporary directory.
}
\seealso{
  \code{\link{tempdir}}
}

