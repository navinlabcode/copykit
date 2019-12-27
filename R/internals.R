#' @title
#' Internal CopyKit functions
#'
#' @description
#' Methods to get or set internal fields from the scCNA class

#' @export
setMethod("segment_ratios", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the segment_ratios data within the assay slot
  SummarizedExperiment::assay(x, "segment_ratios")
})


#' @export
setMethod("ratios", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the ratios data slot
  SummarizedExperiment::assay(x, "ratios")
})

#' @export
setMethod("bin_counts", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the bin_counts data slot
  SummarizedExperiment::assay(x, "bin_counts")
})

#' @export
setMethod("phylo", "scCNA", function(x) {
  # accessor for the phylo slot
  out <- x@phylo
  out
})

#' @export
setReplaceMethod("phylo", "scCNA", function(x, value) {
  # setter method for phylo slot
  x@phylo <- value
  x
})

#' @export
setMethod("distMat", "scCNA", function(x) {
  # accessor for the distMat slot
  out <- x@distMat
  out
})

#' @export
setReplaceMethod("distMat", "scCNA", function(x, value) {
  # setter method for distMat slot
  x@distMat <- value
  x
})

#' @export
setMethod("graph", "scCNA", function(x) {
  # accessor for the getGraph slot
  out <- x@graph
  out
})

#' @export
setReplaceMethod("graph", "scCNA", function(x, value) {
  # setter method for getGraph slot
  x@graph <- value
  x
})

#' @export
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "scCNA", function(object) {
  callNextMethod()
  cat(
    "rowRanges has: ", length(SummarizedExperiment::rowRanges(object)), " ranges\n",
    sep = "",
    "Phylo:"
  )
})
