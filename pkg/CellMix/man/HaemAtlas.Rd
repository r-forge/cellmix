\docType{data}
\name{HaemAtlas}
\alias{HaemAtlas}
\title{HaemAtlas Dataset - Immune Cells}
\source{
  \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2680378/bin/blood-2008-06-162958_TableS5.xls}
}
\description{
  The dataset contains markers for CD4+ and CD8+ T-cells,
  lymphocytes, monocytes, B-cells, NK cells, and
  granulocytes obtained from the HaemAtlas \cite{Watkins et
  al. (2009)}.
}
\details{
  The dataset was created using the internal function
  \code{.HaemAtlas}.

  The probe ids provided in the supplementary data file are
  somehow incorrect, and do not correspond to the Illumina
  ids found in the illuminaHumanv2.db annotation package.
  "Correct" probe ids were recovered from their nucelotide
  sequences, via their nuIDs.
}
\references{
  Watkins Na, Gusnanto A, Bono Bd, De S, Miranda-Saavedra
  D, Hardie DL, Angenent WGJ, Attwood AP, Ellis PD, Erber
  W, Foad NS, Garner SF, Isacke CM, Jolley J, Koch K,
  Macaulay IC, Morley SL, Rendon A, Rice KM, Taylor N,
  Thijssen-Timmer DC, Tijssen MR, Schoot CEvd, Wernisch L,
  Winzer T, Dudbridge F, Buckley CD, Langford CF, Teichmann
  S, Gottgens B and Ouwehand WH (2009). "A HaemAtlas:
  characterizing gene expression in differentiated human
  blood cells." _Blood_, *113*(19), pp. e1-9. ISSN
  1528-0020, <URL:
  http://dx.doi.org/10.1182/blood-2008-06-162958>, <URL:
  http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=2680378\&tool=pmcentrez\&rendertype=abstract>.
}
\seealso{
  Other HaemAtlas: \code{\link{.mHaemAtlas}}
}
\keyword{datasets}

