\name{.qprog}
\alias{.qprog}
\title{Proportion or Signature Estimation by Quadratic Programming}
\usage{
  .qprog(X, seed, data = NULL, exact = TRUE, log = NULL,
    ...)
}
\arguments{
  \item{X}{target matrix}

  \item{seed}{initial NMF model as an \code{\link{NMF}}
  object.}

  \item{data}{markers as \code{\link{MarkerList}} object.}

  \item{exact}{logical that specifies if one should impose
  a sum-to-one (\code{TRUE}) or sum-to-less-than-one
  (\code{FALSE}) constraint on the proportions.}

  \item{log}{indicates if the data should be converted to
  linear scale (not log-trasnformed). If \code{TRUE} the
  data is exponentialised (using base 2). If \code{FALSE}
  the data is left unchanged (the detected log scale is
  displayed in verbose mode). If a number, then it is used
  as the base to exponentialise the data.}

  \item{...}{extra parameters passed to
  \code{\link{mlsei}}.}
}
\value{
  an \code{\link{NMF}} object.
}
\description{
  This is the workhorse function for the ged algorithms
  \sQuote{qprog}, \sQuote{cs-qprog}, and \sQuote{DSA}.
}
\details{
  The function implements deconvolution methods that use
  quadratic programming techniques to estimate proportions
  from signatures (algorithm
  \code{\link[=qprog-ged]{qprog}} from \cite{Gong et al.
  (2011)}), or both signatures and proportions from marker
  genes (algorithm \code{\link[=DSA-ged]{DSA}} from
  \cite{Zhong et al. (2013)}).
}
\references{
  Gong T, Hartmann N, Kohane IS, Brinkmann V, Staedtler F,
  Letzkus M, Bongiovanni S and Szustakowski JD (2011).
  "Optimal deconvolution of transcriptional profiling data
  using quadratic programming with application to complex
  clinical blood samples." _PloS one_, *6*(11), pp. e27156.
  ISSN 1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0027156>, <URL:
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=3217948\&tool=pmcentrez\&rendertype=abstract>.

  Zhong Y, Wan Y, Pang K, Chow LM and Liu Z (2013).
  "Digital sorting of complex tissues for cell
  type-specific gene expression profiles." _BMC
  Bioinformatics_, *14*(1), pp. 89. ISSN 1471-2105, <URL:
  http://dx.doi.org/10.1186/1471-2105-14-89>, <URL:
  http://www.biomedcentral.com/1471-2105/14/89>.
}
\keyword{internal}

