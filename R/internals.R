#' @title
#' Internal CopyKit functions
#'
#' @description
#' Methods to get or set internal fields from the scCNA class

#' export
setMethod("segment_ratios", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the segment_ratios data within the assay slot
  SummarizedExperiment::assay(x, "segment_ratios")
})
