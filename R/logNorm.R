#' logNorm()
#'
#' Computes a log transformation of the selected assay
#'
#' @param scCNA scCNA object.
#' @param transform String specifying the transformation to apply to the selected
#' assay.
#' @param assay String with the name of the assay to pull data from to run the
#' segmentation.
#' @param name String with the name for the target slot for the resulting
#' transformed counts.
#'
#' @return A data frame with log transformed counts inside the
#' \code{\link[SummarizedExperiment]{assay}} slot.
#'
#' @importFrom SummarizedExperiment assay
#'
#' @export
#'
#' @examples
logNorm <- function(scCNA,
                    transform = c("log", "log2", "log10", "log1p"),
                    assay = 'segment_ratios',
                    name = 'logr') {

  transform <- match.arg(transform)

  # obtaining data
  seg_ratios <- SummarizedExperiment::assay(scCNA, assay)

  #saving logr
  seg_ratios[seg_ratios == 0] <- 1e-3

  if (transform == 'log') {
    seg_ratios_logr <- log(seg_ratios)
  } else if (transform == 'log2') {
    seg_ratios_logr <- log2(seg_ratios)
  } else if (transform == 'log1p') {
    seg_ratios_logr <- log1p(seg_ratios)
  } else if (transform == 'log10') {
    seg_ratios_logr <- log10(seg_ratios)
  }

  SummarizedExperiment::assay(scCNA, name) <- round(seg_ratios_logr, 2)

  return(scCNA)

}
