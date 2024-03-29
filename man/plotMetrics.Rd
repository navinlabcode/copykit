% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMetrics.R
\name{plotMetrics}
\alias{plotMetrics}
\title{plotMetrics}
\usage{
plotMetrics(scCNA, metric, label = NULL, dodge.width = 0, ncol = 2)
}
\arguments{
\item{scCNA}{scCNA object.}

\item{metric}{A character indicating which elements of \code{colData()}
should be plotted.}

\item{label}{A character indicating which element of the \code{colData()} to
color the plots.}

\item{dodge.width}{A numeric that adds dodge between the label categories.}

\item{ncol}{A Integer specifying the number of columns to be used for the
panels of a multi-facet plot.}
}
\value{
A ggplot object with swarm plots of the selected metrics.
}
\description{
Plots swarm plots from selected elements of \code{\link[SummarizedExperiment]{colData}}.
}
\examples{
copykit_obj <- mock_bincounts()
copykit_obj <- runMetrics(copykit_obj)
}
\author{
Darlan Conterno Minussi
}
