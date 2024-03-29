% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotScree.R
\name{plotScree}
\alias{plotScree}
\title{plotScree}
\usage{
plotScree(scCNA, ncomponents = 20)
}
\arguments{
\item{scCNA}{The CopyKit object}

\item{ncomponents}{Number of principal components to plot.}
}
\value{
A ggplot object with The variance explained per principal component.
}
\description{
Plots the variance explained by the different principal components
}
\examples{
set.seed(1000)
copykit_obj <- copykit_example_filtered()[,sample(400)]
copykit_obj <- runPca(copykit_obj)
plotScree(copykit_obj)

}
