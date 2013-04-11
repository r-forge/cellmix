\name{xsanitize}
\alias{xsanitize}
\title{Partial Sanitizer for xtable}
\usage{
  xsanitize(skip)
}
\arguments{
  \item{skip}{character string that contain the character
  not to sanitize.}
}
\value{
  a function that can be used in argument
  \code{sanitize.text.function} of \code{print.xtable}.
}
\description{
  Text sanitizing function to use in argument
  \code{sanitize.text.function} of
  \code{\link{print.xtable}}. This function is used in
  \pkg{CellMix} to format the content internal registries
  for inclusion in vignettes, as they contain citation
  commands (\\cite{}).
}
\examples{
library(xtable)
print(xtable(gedAlgorithmInfo(show=FALSE), citet=TRUE), sanitize.text.function=xsanitize("\\\\{}$"))
}
\seealso{
  \code{\link{print.xtable}}
}
\keyword{internal}

