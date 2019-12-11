#' Creates UMAP embedding
#'
#' Creates a umap embedding using the package uwot from the segment ratios values
#'
#' @param scCNA scCNA object.
#' @param seed Sets a seed for the pseudorandom number generator.
#' @param ... Additional parameters passed to \code{uwot::umap}.
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object.
#' @export
#'
#' @examples

runUmap <- function(scCNA,
                    seed = 17,
                    ...) {
  seg_data <- t(segment_ratios(scCNA))

  set.seed(seed)

  dat_umap <- uwot::umap(seg_data,
                         ...)

  SingleCellExperiment::reducedDims(scCNA) <- list(umap = dat_umap)

  return(scCNA)

}
