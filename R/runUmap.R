#' Creates UMAP embedding
#'
#' Creates a umap embedding using the package uwot from the segment ratios
#'  values
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param assay String with the name of the assay to pull data from to make the
#' embedding.
#' @param seed Sets a seed for the pseudorandom number generator.
#' @param name String specifying the name to be used to store the result in the
#' reducedDims of the output.
#' @param min_dist The effective minimum distance between embedded points.
#' Smaller values will result in a more clustered/clumped embedding where nearby
#' points on the manifold are drawn closer together, while larger values will
#' result on a more even dispersal of points. The value should be set relative
#' to the spread value, which determines the scale at which embedded points
#' will be spread out. See \code{\link[uwot]{umap}}.
#' @param n_neighbors The size of local neighborhood (in terms of number of
#' neighboring sample points) used for manifold approximation.
#' Larger values result in more global views of the manifold,
#' while smaller values result in more local data being preserved.
#' In general values should be in the range 2 to 100.
#' See \code{\link[uwot]{umap}}.
#' @param ncomponents The dimension of the space to embed into. See
#' \code{\link[uwot]{umap}}.
#' @param ... Additional parameters passed to \code{\link[uwot]{umap}}.
#'
#' @importFrom uwot umap
#' @importFrom SummarizedExperiment assay
#' @importFrom withr with_seed
#'
#' @return A reduced dimension representation with UMAP in the slot
#'  \code{reducedDim} from scCNA object. Access reduced dimensions slot with:
#'  \code{reducedDim(scCNA, 'umap', withDimnames = FALSE)}
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- runUmap(copykit_obj)
runUmap <- function(scCNA,
                    assay = "logr",
                    seed = 17,
                    min_dist = 0,
                    n_neighbors = 50,
                    name = "umap",
                    ncomponents = 2,
                    ...) {
    seg_data <- t(SummarizedExperiment::assay(scCNA, assay)) %>%
        as.data.frame()

    message("Using assay: ", assay)
    message("Embedding data with UMAP. Using seed ", seed)
    withr::with_seed(
        seed,
        dat_umap <- uwot::umap(seg_data,
                               min_dist = min_dist,
                               n_neighbors = n_neighbors,
                               ncomponents = n_components,
                               ...)
    )

    SingleCellExperiment::reducedDim(scCNA, type = name) <- dat_umap

    message(
        "Access reduced dimensions slot with: reducedDim(scCNA, 'umap')."
    )
    message("Done.")

    return(scCNA)
}
