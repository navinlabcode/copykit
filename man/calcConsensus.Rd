% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcConsensus.R
\name{calcConsensus}
\alias{calcConsensus}
\title{Calculate a consensus matrix of segment means based on \code{colData}}
\usage{
calcConsensus(
  scCNA,
  assay = "segment_ratios",
  consensus_by = "subclones",
  fun = c("median", "mean"),
  BPPARAM = bpparam()
)
}
\arguments{
\item{scCNA}{The CopyKit object.}

\item{assay}{String with the name of the assay to pull data from to calculate
the consensus matrix.}

\item{consensus_by}{A string with the column from colData that will be used
to isolate the cells by factor and calculate the consensus.}

\item{fun}{A string indicating the summarizing function to be used.}

\item{BPPARAM}{A \linkS4class{BiocParallelParam} specifying how the function
should be parallelized.}
}
\value{
A consensus matrix stored in the consensus slot of the CopyKit object
}
\description{
Calculate a consensus matrix of segment means based on \code{colData}
}
\details{
Consensus profiles are calculated by averaging or taking the median
 of the ith segment mean of all single cells assigned to the same element of
 \link{colData},
}
\examples{
copykit_obj <- copykit_example_filtered()
copykit_obj <- findClusters(copykit_obj)
copykit_obj <- calcConsensus(copykit_obj)
}
