#' Creates UMAP embedding
#'
#' Creates a umap embedding using the package uwot from the segment ratios values
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param assay String with the name of the assay to pull data from to make the embedding.
#' @param seed Sets a seed for the pseudorandom number generator.
#' @param name String specifying the name to be used to store the result in the
#' reducedDims of the output.
#' @param ... Additional parameters passed to \code{uwot::umap}.
#'
#' @importFrom uwot umap
#' @importFrom SummarizedExperiment assay
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object. Access reduced dimensions slot with: \code{SingleCellExperiment::reducedDim(scCNA, 'umap', withDimnames = FALSE)}
#' @export
#'
#' @examples

runUmap <- function(scCNA,
                    assay = "segment_ratios",
                    seed = 17,
                    name = "umap",
                    ...) {
  seg_data <- t(SummarizedExperiment::assay(scCNA, assay)) %>%
    as.data.frame()

  message(paste("Using assay:", assay))
  message(paste("Embedding data with UMAP. Using seed", seed))
  set.seed(seed)

  dat_umap <- uwot::umap(seg_data,
                         ...)

  SingleCellExperiment::reducedDim(scCNA, type = name) <- dat_umap

  message(
    "Access reduced dimensions slot with: SingleCellExperiment::reducedDim(scCNA, 'umap')."
  )
  message("Done.")

  return(scCNA)

}
