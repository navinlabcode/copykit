---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
# CopyKit

<!-- badges: start -->
<!-- badges: end -->

  **CopyKit** provides a toolkit for the analysis of single-cell copy number datasets. It includes functions to read data from [Varbin](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4417119/) or [10X CNA datasets](https://www.10xgenomics.com/solutions/single-cell-cnv/).
  
  A common workflow with **CopyKit** consists in reading the dataset and using to filter the noisy cells out, clustering and plotting heatmaps for further analysis.


## Installation

You can install the development version of CopyKit from github with:

``` r
devtools::install_github("navinlabcode/copykit")
```

## Tutorial

