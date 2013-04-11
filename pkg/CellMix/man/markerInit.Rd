\name{markerInit}
\alias{markerInit}
\title{Initialisation Method for semi-NMF Algorithms}
\usage{
  markerInit(method, v, x, data = NULL,
    markers = c("prior+semi", "semi"), ratio = NULL,
    eps = NULL, ...)
}
\arguments{
  \item{method}{method name, used in case of error only}

  \item{v}{target matrix}

  \item{x}{initial NMF model}

  \item{data}{marker list}

  \item{ratio}{expression ratio of markers between its cell
  type and other cell types.}

  \item{eps}{lower limit for cell type signature expression
  levels or proportions. This is used to avoid null values
  that might be problematic for some deconvolution
  algorithm, e.g., for the method \sQuote{ssKL} that is
  based on the KL-divergence, which computes logs.}

  \item{markers}{indicates what the markers are used for:
  \describe{ \item{sQuote{prior}}{uses
  \code{\link[=DSA-ged]{DSA}} proportion estimation method
  from \cite{Zhong2013} to compute sensible initial
  proportions from average marker expression profiles in
  the mixed sample data.} \item{sQuote{semi}}{enforces
  marker block patterns after each iteration.}
  \item{sQuote{post}}{\emph{a posteriori} assigns estimated
  signatures;} }}

  \item{...}{extra parameters currently not used}
}
\value{
  the modified NMF model.
}
\description{
  Pre-process and attach markers to an NMF model so that
  they can be enforced at each iteration. It adds a static
  variable 'markerIdx' to the workspace.
}
\keyword{internal}

