#' Creates UMAP embedding
#'
#' Creates a umap embedding using the package uwot from the segment ratios values
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param seed Sets a seed for the pseudorandom number generator.
#' @param ... Additional parameters passed to \code{uwot::umap}.
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object. Access reduced dimensions slot with: \code{SingleCellExperiment::reducedDims(scCNA, 'umap')}
#' @export
#'
#' @examples

runUmap <- function(scCNA,
                    seed = 17,
                    ...) {
  seg_data <- t(segment_ratios(scCNA))

  message(paste("Embedding data with UMAP. Using seed:", seed))
  set.seed(seed)

  dat_umap <- uwot::umap(seg_data,
                         ...)

  SingleCellExperiment::reducedDims(scCNA) <- list(umap = dat_umap)

  message("Done. Access reduced dimensions slot with: SingleCellExperiment::reducedDims(scCNA, 'umap')")

  return(scCNA)

}
