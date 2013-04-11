\name{markerScoreAbbas}
\alias{markerScoreAbbas}
\title{Marker Scoring Method: Abbas et al. (2009)}
\usage{
  markerScoreAbbas(object, data, statistic = "p.value",
    ntop = 2, log = !is_logscale(object), lbase = 2,
    vsall = FALSE)
}
\arguments{
  \item{statistic}{statistic to use as a score. The method
  computes the following quantities for each comparison
  between the most expressing cell type (highest
  within-group mean) with the i-th most expressing other
  cell type: \describe{ \item{p.value<i>}{p-value
  (t-test).} \item{dm<i>}{difference in means.}
  \item{statistic<i>}{statistic value (t statistic).}
  \item{dmM<i>}{minimum difference = min(top) -
  max(top<i>).} \item{fold<i>}{fold change in means.}
  \item{mMfold<i>}{minimum fold change =
  min(top)/max(top<i>).} \item{top<i>}{group index of the
  i-th most expressing cell type. It is relative to the
  levels of \code{data}.} }}

  \item{ntop}{the number of groups for which statistics
  should be computed (>= 2). If \code{ntop=2}, only
  statistics between the highest and second-highest
  expressing cell-type are computed. If \code{ntop=3},
  statistics between the highest and third-highest are also
  computed, and so on. Use \code{ntop=Inf} to compare the
  highest-expressing group to all other groups.}

  \item{log}{a logical specifying if the should the data be
  log-transformed before computing the p-values. Default is
  to log-transform the data if is not already in log-scale,
  which is determined by the function
  \code{\link{is_logscale}}.}

  \item{lbase}{log base to used if the data is
  log-transformed.}

  \item{vsall}{logical used when \code{ntop=2} that
  indicates if the comparison should be carried out between
  the highest-expressing cell-type and the rest of all
  other cell-types.}

  \item{object}{A numeric \code{matrix}, an object of class
  \code{\link[Biobase]{ExpressionSet}}, or a
  \code{\linkS4class{MarkerList}} object.}

  \item{data}{main extra data used by the scoring and
  selection methods.

  If \var{object} is a matrix-like object, \var{data} is
  generally a factor or grouping variable, as a vector,
  that defines the cell-type for each -- pure -- sample. If
  a vector, this argument is converted into a factor with
  the levels in their order of appearance in \code{data},
  by \code{factor(data, levels=unique(data))}. This is to
  obtain levels in an order that is consistent with the
  samples' order.

  If \var{object} is a \code{MarkerList} object, then
  \var{data} is generally a matrix-like object that
  contains expression data.}
}
\value{
  \item{extractMarkers}{a numeric matrix with the following
  named columns: \itemize{ \item{top}{the index of the
  highest-expressing cell-type, as defined by the levels of
  the \code{factor} derived from argument \code{group}.}
  \item{p.value2}{ p-value from the comparison between the
  highest-expressing (top) and second highest-expressing
  (second top) cell-type.} \item{dm2}{difference in mean
  expression between the top and second top cell-type}
  \item{dmM2}{minimum difference of expression between the
  top and second top cell-type} \item{fold2}{fold change
  mean expression between the top and second top cell-type}
  \item{mMfold2}{minimum fold change expression between the
  top and second top cell-type} \item{dm3, mdM3, fold3,
  mMfold3}{ same quantitities computes between the top and
  third-top cell-type.} }

  The result matrix has an attribute \code{'types'} that
  contains the levels of the original (or converted) factor
  \code{group}.}
}
\description{
  \code{markerScoreAbbas} implements the scoring/selection
  method proposed by \cite{Abbas et al. (2009)}, to select
  marker genes from pure cell type samples.
}
\details{
  The method \sQuote{Abbas} uses a t-test approach, so that
  the data is assumed to contain at least 2 pure samples
  per cell-type. It implements the method from \cite{Abbas
  et al. (2009)}:

  "[...] top differentially expressed (based on 95% fold
  change confidence intervals from Student's T-test)
  probesets were determined by comparing each probe's
  highest-expressed group with the next highest-expressed
  group in order to find probesets that are good markers
  for each cell population. This step was repeated with
  comparison between the top group and the third-highest
  group in order to also include probesets that were strong
  markers for two cell populations."

  For each gene, the highest-expressing cell type is
  determined by ordering them by mean expression.
  Comparisons and p-value computations are performed using
  the fast t-test implementation from
  \code{\link{rowttests}} in the \pkg{genefilter} package.
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
  Other markerScore: \code{\link{extractMarkers}},
  \code{\link{markerScoreHSD}},
  \code{\link{markerScoreMaxcol}},
  \code{\link{markerScoreMethod}},
  \code{\link{markerScoreScorem}},
  \code{\link{scoreMarkers}},
  \code{\link{selectMarkers.markerScore_HSD}},
  \code{\link{selectMarkers.markerScore_scorem}}
}

