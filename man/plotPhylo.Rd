% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotPhylo.R
\name{plotPhylo}
\alias{plotPhylo}
\title{plotPhylo()}
\usage{
plotPhylo(
  scCNA,
  label = NULL,
  label_colors = NULL,
  consensus = FALSE,
  group = NULL
)
}
\arguments{
\item{scCNA}{scCNA object.}

\item{label}{A string with the element of
\code{\link[SummarizedExperiment]{colData}}. to annotate the tips of the
tree.}

\item{label_colors}{A named list with colors for the label annotation.
Must match label length}

\item{consensus}{A boolean indicating if the consensus phylogenetic tree
should be plotted.}

\item{group}{A string that if provided will plot the tip labels as pie charts
with the proportions from the provided element from
\code{\link[SummarizedExperiment]{colData}}}
}
\value{
A ggplot object with a phylogenetic tree visualization.
}
\description{
Plots a phylogenetic tree.
}
\examples{
set.seed(1000)
copykit_obj <- copykit_example_filtered()[,sample(100)]
copykit_obj <- findClusters(copykit_obj)
copykit_obj <- runPhylo(copykit_obj)
plotPhylo(copykit_obj, label = "subclones")
}
\author{
Junke Wang
}
