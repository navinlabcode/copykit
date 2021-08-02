#' Run metrics
#'
#' Calculates the overdispersion and the breakpoint counts for each cell.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details Adds the metrics to the scCNA \code{\link[SummarizedExperiment]{colData}}.
#'  Those metrics can be used for subsetting the data if desired.
#'  results can be visualized with \code{\link{plotMetrics}}.
#'
#' @return Adds columns 'overdispersion' and 'breakpoint_count' to
#'  \code{\link[SummarizedExperiment]{colData}}.
#'
#' @export
#' @import ggplot2
#'
#' @examples
#'
#'
runMetrics <- function(scCNA,
                       BPPARAM = bpparam()) {

  ###################
  # Retrieving data

  dat_seg <- segment_ratios(scCNA)
  dat_rat <- ratios(scCNA)
  dat_bin <- bin_counts(scCNA)
  rg <- SummarizedExperiment::rowRanges(scCNA)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:01 2021
  # overdispersion ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:10 2021

  message("Calculating overdispersion.")

  overdisp <- BiocParallel::bplapply(dat_bin,
                                  overdispersion,
                                  BPPARAM = BPPARAM)

  overdisp <- unlist(overdisp)

  SummarizedExperiment::colData(scCNA)$overdispersion <- overdisp

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:25 2021
  # Breakpoint count
  # Performed for every chromosome
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:34 2021

  scCNA <- .countBreakpoints(scCNA)

  message("Done.")

  return(scCNA)

}

