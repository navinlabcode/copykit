#' Run phylogenetic analysis
#'
#' Performs phylogenetic analysis using the segment ratios data
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param method Phylogenetic method to be run, currently only accepts "nj" (neighbor-joining).
#' @param metric distance metric passed to construct the phylogeny.
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `amap::Dist`
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object. Access reduced dimensions slot with: \code{SingleCellExperiment::reducedDim(scCNA, 'umap', withDimnames = FALSE)}
#' @export
#'
#' @examples

runPhylo <- function(scCNA,
                     method = "nj",
                     metric = "euclidean",
                     n_threads = 1) {

  seg_data <- t(segment_ratios(scCNA)) %>%
    as.data.frame()

  if (nrow(seg_data) > 500) {
    message(paste("Your dataset has:",
                  nrow(seg_data), "Cells"))
    message("Calculating the distance matrix may take a long time with large number of cells")
    message("Set number of threads with n_threads for parallel processing if possible to speed up.")
  }

  # ordering cells
  tree <- ape::nj(amap::Dist(seg_data, nbproc = n_threads))

  tree <- ape::ladderize(tree)

  phylo(scCNA) <- tree

  message("Access slot with copykit::phylo(scCNA).")
  message("Done.")
  return(scCNA)

}

