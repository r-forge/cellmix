\docType{methods}
\name{gedBlood}
\alias{asCBC}
\alias{asCBC,character-method}
\alias{asCBC,MarkerList-method}
\alias{asCBC-methods}
\alias{asCBC,NMF-method}
\alias{gCBC}
\alias{gedBlood}
\alias{refCBC}
\title{Blood Sample Deconvolution}
\format{Named num [1:5] 0.005 0.295 0.03 0.57 0.1
 - attr(*, "names")= chr [1:5] "Basophils" "Lymphocytes" "Eosinophils" "Neutrophils" ...}
\usage{
  gedBlood(object, ..., normalize = TRUE, verbose = FALSE)

  asCBC(object, ...)

  \S4method{asCBC}{character}(object, drop = FALSE,
    quiet = FALSE)

  \S4method{asCBC}{NMF}(object, drop = TRUE, ...)

  refCBC

  gCBC(n = 1, sampleNames = NULL, counts = NULL)
}
\arguments{
  \item{object}{target data, specified in any format
  supported by \code{\link{ged}}. For \code{asCBC}, an
  object for with suitable \code{asCBC} method defined.}

  \item{...}{extra arguments passed to
  \code{\link{gedProportions}}.}

  \item{drop}{logical that indicates if elements in
  \code{object} that cannot be mapped to a cell type should
  be removed from the returned mapping.}

  \item{quiet}{logical that indicates that the mapping
  should be performed quietly. If \code{FALSE}, then an
  error is thrown if none of the elements can be mapped,
  or, if in addition \code{drop=FALSE}, a warning is thrown
  if only some of the elements could be mapped.}

  \item{n}{number of samples in the generated CBC matrix}

  \item{sampleNames}{names of the samples, recycled or
  truncated if necessary, to match \code{n}.}

  \item{counts}{CBC data to use instead of the defaults. It
  must be a numeric vector.}

  \item{normalize}{character string that specifies the
  normalisation method to apply jointly to the combined
  data (signatures,data). The normalisation is performed
  after transforming the data and/or signatures if
  necessary.}

  \item{verbose}{logical that toggles verbosity. A number
  (integer) can be passed to specify the verbosity level
  (the higher the more messages are output). Passing
  \code{verbose=Inf} toggles debug mode (all messages).
  Note that because it appears after \code{...} it must be
  fully named.}
}
\description{
  The functions described here are dedicated to gene
  expression deconvolution of blood samples (i.e. whole
  blood or PBMCs).

  \code{gedBlood} uses the methodology defined by
  \cite{Abbas et al. (2009)}, which uses a fixed set of 17
  cell type-specific signatures to estimate cell
  proportions in blood samples. Each signature corresponds
  to a white blood cell in resting or activated state (See
  section Details).

  \code{asCBC} has methods defined for \code{NMF} models
  and \code{Markerlist} objects. See each method's
  description for more details.

  \code{refCBC} is a numeric vector that contains average
  Complete Blood Count proportions (CBC) in healthy
  persons, based on empirical studies in healthy patients.
  It contains proportions for Basophils, Lymphocytes,
  Eosinophils, Neutrophils and Monocytes.

  \code{gCBC} generates a matrix of average Complete Blood
  Count proportions (CBC) for a given number of samples.
  The default proportions are based on empirical studies in
  healthy patients (see \code{\link{refCBC}}), and each
  sample get assigned the same proportions.
}
\details{
  The signatures used by \code{gedBlood} were designed by
  \cite{Abbas et al. (2009)} to optimise their
  deconvolution power. They are available in the
  \pkg{CellMix} as dataset \code{\link{Abbas}}.

  \code{gedBlood} is essentially a shortcut for
  \code{gedProportions(object, Abbas, ...)}, see
  \code{\link{gedProportions}} for details on the possible
  arguments.

  Currently \code{asCBC} methods will correctly work only
  on objects that have cell types that match exactly names
  of signatures in the \code{\link{Abbas}} dataset.
}
\section{Methods}{
  \describe{

  \item{asCBC}{\code{signature(object = "character")}: This
  is the workhorse method that maps immune/blood cell type
  names to the CBC cell types: Monocytes, Basophils,
  Lymphocytes, Neutrophils and Eosinophils.

  It returns a factor, whose names are elements of
  \code{object} and the values are their corresponding CBC
  cell type. If \code{drop=FALSE} the result is of the same
  length as \code{object}, otherwise it only contains
  elements that could be mappped to a cell type. }

  \item{asCBC}{\code{signature(object = "NMF")}: The result
  of gene expression deconvolution performed by
  \code{\link{ged}} are stored in
  \code{\linkS4class{NMFstd}} model objects, which contain
  the cell type-specific signatures and/or cell relative
  proportions.

  This method aggregates, i.e. sums up, the cell
  proportions and averages the signatures of cell types
  from each of the CBC groups that are available in the
  data. }

  \item{asCBC}{\code{signature(object = "MarkerList")}:
  This method combines markers of cell types that belong to
  the same CBC group. }

  }
}
\examples{
# load kidney transplant data (Shen-Orr et al. (2010))
eset <- ExpressionMix('GSE20300')

# true CBC proportions are stored as -- mixture -- coefficients
cbc <- coef(eset)
head(cbc)
profplot(cbc, lab='')

# devonvolve using basis signature matrix from Abbas et al. (2009)
res <- gedBlood(eset, verbose=TRUE)
# estimated proportions are stored as -- mixture -- coefficients
p <- coef(res)
str(p)
# the Abbas basis matrix includes detailed immune cell types
rownames(p)
# or: basisnames(res)

# aggregate into CBC
ecbc <- asCBC(res)
# plot Estimated vs. Measured CBC
profplot(cbc, ecbc)
# default proportions
gCBC()

# for 5 samples
gCBC(5)

# setting sample names
gCBC(5, sampleNames = letters[1:10]) # names are truncated if necessary
}
\references{
  Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z and Clark
  HF (2009). "Deconvolution of blood microarray data
  identifies cellular activation patterns in systemic lupus
  erythematosus." _PloS one_, *4*(7), pp. e6098. ISSN
  1932-6203, <URL:
  http://dx.doi.org/10.1371/journal.pone.0006098>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/19568420>.
}
\seealso{
  \code{\link{gedProportions}}, \code{\link{ged}}
}
\keyword{datasets}
\keyword{methods}

