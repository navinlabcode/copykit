% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runPhylo.R
\name{runPhylo}
\alias{runPhylo}
\title{Run phylogenetic analysis}
\usage{
runPhylo(
  scCNA,
  method = "nj",
  metric = "euclidean",
  assay = "segment_ratios",
  n_threads = parallel::detectCores()/4
)
}
\arguments{
\item{scCNA}{scCNA object.}

\item{method}{Phylogenetic method to be run, currently accepts "nj" (neighbor-joining) and "me" (minimum evolution). Defaults to "nj".}

\item{metric}{distance metric passed to construct the phylogeny (Defaults to "euclidean").}

\item{assay}{String with the name of the assay to pull data from to run phylogenetic analysis. Note that only assay named "integer" will be treated as integer.}

\item{n_threads}{Number of threads used to calculate the distance matrix. Passed to `amap::Dist`}
}
\value{
A rooted phylogenetic tree object in the slot \code{phylo} from scCNA object. Access phylo slot with: \code{copykit::phylo(scCNA)}
}
\description{
Performs phylogenetic analysis
}
\examples{
set.seed(1000)
copykit_obj <- copykit_example_filtered()[,sample(50)]
copykit_obj <- runPhylo(copykit_obj)
}
\author{
Darlan Conterno Minussi

Junke Wang
}
