\name{nuIDdecode}
\alias{nuIDdecode}
\title{Convert nuID to Nucleotide Sequence}
\source{
  function \code{lumi::id2seq} in the \pkg{lumi} package
}
\usage{
  nuIDdecode(id, error = TRUE)
}
\arguments{
  \item{id}{nuID character vector}

  \item{error}{a logical indicating whether an error should
  be thrown in case of invalid input nuID (default) or only
  a warning, as in the original function
  \code{lumi::id2seq}.}
}
\description{
  The function \code{nuIDdecode} converts a nuID string
  back to the nucleotide sequence it encodes. nuIDs are
  identifiers used as primary keys in some Illumina
  annotation packages, and are based on a hash of the probe
  sequence itself.
}
\details{
  This function is an adapted version of the
  \code{lumi::id2seq} from the \pkg{lumi} package, so that
  it can throws errors instead of warnings. It is used in
  \code{\link{idtype}} to infer the type of nuID vectors.
}
\examples{
nuIDdecode('XERDqYYc2A')
try(nuIDdecode('XERDqYYc2F'))
nuIDdecode('XERDqYYc2F', error=FALSE)
}
\author{
  Pan Du

  Adaptation for CellMix: Renaud Gaujoux
}
\references{
  Du P, Kibbe Wa and Lin SM (2007). "nuID: a universal
  naming scheme of oligonucleotides for illumina,
  affymetrix, and other microarrays." _Biology direct_,
  *2*, pp. 16. ISSN 1745-6150, <URL:
  http://dx.doi.org/10.1186/1745-6150-2-16>, <URL:
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=1891274\&tool=pmcentrez\&rendertype=abstract>.
}

