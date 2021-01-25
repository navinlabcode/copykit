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
#' @importFrom uwot umap
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object. Access reduced dimensions slot with: \code{SingleCellExperiment::reducedDim(scCNA, 'umap', withDimnames = FALSE)}
#' @export
#'
#' @examples

runUmap <- function(scCNA,
                    seed = 17,
                    ...) {
  seg_data <- t(segment_ratios(scCNA)) %>%
    as.data.frame()

  message(paste("Embedding data with UMAP. Using seed", seed))
  set.seed(seed)

  dat_umap <- uwot::umap(seg_data,
                         ...)
  SingleCellExperiment::reducedDims(scCNA) <- list(umap = dat_umap)

  message(
    "Access reduced dimensions slot with: SingleCellExperiment::reducedDim(scCNA, 'umap')."
  )
  message("Done.")

  return(scCNA)

}
