% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sf.R
\name{computeSizeFactors}
\alias{computeSizeFactors}
\title{computeSizeFactors}
\usage{
computeSizeFactors(ggd, factorGroups = NULL)
}
\arguments{
\item{ggd}{A GenoGAMDataSet object.}

\item{factorGroups}{A list of grouped IDs (same as the colnames of
the GenoGAMDataSet object). Each element of the list represents
a group of samples within which size factors are computed.
If NULL all samples are regarded to belong to one group.}
}
\value{
A GenoGAMDataSet object, where the sizeFactors slot is updated.
}
\description{
The function computes the size factors for given factor groups based on
the DESeq2 package.
}
\examples{
ggd <- makeTestGenoGAMDataSet()
ggd <- computeSizeFactors(ggd)
sizeFactors(ggd)
groups <- list(c("wt_1", "wt_2"), c("mutant_1", "mutant_2"))
ggd <- computeSizeFactors(ggd, factorGroups = groups)
sizeFactors(ggd)
}
\author{
Georg Stricker \email{georg.stricker@in.tum.de}
}
