\name{DSection}
\alias{DSection}
\title{DSection Gene Expression Deconvolution Method}
\source{
  \url{http://www.cs.tut.fi/~erkkila2/software/dsection/DSection.m}
}
\usage{
  DSection(Y, p0, nSamples, groups = NULL,
    nBurnIn = 4 * nSamples, W0 = max(nrow(p0), 10),
    W_proposal = 100, samplep = TRUE, summarize = TRUE,
    verbose = FALSE)
}
\arguments{
  \item{Y}{<I-by-J> matrix of measurements from
  heterogeneous tissues. I is the number of
  probes/genes/etc., and J is the number of tissues.}

  \item{p0}{<T-by-J> matrix of prior predictions on cell
  type proportions. T is the number of cell types, and
  columns in \code{p0} must be positive and add up to one.}

  \item{groups}{<1-by-J> vector of treatment indices, so
  that unique(Treatment) = [1,2,...,C], where C is the
  number of treatments including control, i.e., "no
  treatment", if available.}

  \item{W0}{Prior prediction weight, i.e., degree of
  confidence, on \code{p0}. Defines the peakedness of
  Dirichlet density around p0. NOTE: keep W0 >= T.}

  \item{W_proposal}{Transition kernel weight, defines the
  peakedness of Dirichlet density around p*, the old value.
  The higher W_proposal is, the smaller the proposal steps
  around p* are.}

  \item{nBurnIn}{Amount of burn-in. NOTE: keep nBurnIn >
  0.}

  \item{nSamples}{Amount of sampling. NOTE: keep nSamples >
  0.}

  \item{samplep}{logical value, indicating whether to
  sample from the posterior for cell type proportions
  (\code{TRUE}) or not (\code{FALSE}). SUGGESTED USE:
  sample from the posterior (samplep = 1).}

  \item{summarize}{logical indicating whether only average
  values should be returned -- and computed.}

  \item{verbose}{logical that indicates if verbose messages
  should be shown.}
}
\value{
  A list with the following elements:

  \item{MCData}{ results from the MCMC estimation.}
  \item{x_LS}{ Standard-least square estimate.}
  \item{groups}{ factor defining the groups of samples, if
  any was provided.} \item{call}{ the call to DSection.}
  \item{parameters}{ a list of some of the parameters used
  in the estimation.} \item{p0}{ initial prior on
  proportions.}
}
\description{
  The \emph{DSection} algorithm performs gene expression
  deconvolution when priors on proportions are available,
  using a Markov Chain Monte Carlo approach \cite{Erkkila
  et al. (2010)}.
}
\details{
  In \pkg{CellMix}, this method is registered with the key
  \code{'DSection'}, and is can be applied to gene
  expression data via the function \code{\link{ged}}.

  This function uses the \pkg{RcppOctave} package to run
  the original Matlab code in \emph{Octave}. The
  documentation was extracted from the Matlab source file,
  that can be found in the \pkg{CellMix} package
  "scripts/DSection" subdirectory.

  The Matlab code requires the Octave packages
  \emph{statistics} and \emph{optim} to run properly. These
  packages can be downloaded from Octave-forge:

  \url{http://sourceforge.net/projects/octave/files/Octave\%20Forge\%20Packages/Individual\%20Package\%20Releases/}

  and installed in Octave with:

  \code{pkg install '<path/to/package/tar/gz/file>'}

  or in R

  \code{install.extras('CellMix', 'octave:DSection')}
}
\examples{
#' # random global expression
x <- rmix(3, 50, 10)
dim(x)

# extract true proportions
p <- coef(x)
# add noise to proportions
p0 <- scoef(abs(p + rmatrix(p, dist=rnorm, sd=0.15)))
# check how noisy this makes the proportion estimates
profplot(p, p0)

\dontrun{
# Requires some extra Octave packages to be installed (it defines `gamfit`)
res <- DSection(exprs(x), p0, nSamples=100, nBurnIn=1L)
profplot(p, res$MCData$p)
}
}
\author{
  Original Matlab code: Timo Erkkila

  Wrapper function: Renaud Gaujoux
}
\references{
  \url{http://informatics.systemsbiology.net/DSection/}
  \url{http://www.cs.tut.fi/~erkkila2/software/dsection/index.html}

  Erkkila T, Lehmusvaara S, Ruusuvuori P, Visakorpi T,
  Shmulevich I and Lahdesmaki H (2010). "Probabilistic
  analysis of gene expression measurements from
  heterogeneous tissues." _Bioinformatics (Oxford,
  England)_, *26*(20), pp. 2571-7. ISSN 1367-4811, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btq406>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/20631160
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2951082\&tool=pmcentrez\&rendertype=abstract>.
}

