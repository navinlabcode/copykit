#' Run phylogenetic analysis
#'
#' Performs phylogenetic analysis using the segment ratios data
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param method Phylogenetic method to be run, currently only accepts "nj" (neighbor-joining).
#' @param metric distance metric passed to construct the phylogeny (Defaults to "euclidean").
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `amap::Dist`
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object. Access reduced dimensions slot with: \code{SingleCellExperiment::reducedDim(scCNA, 'umap')}
#' @export
#'
#' @importFrom ape nj ladderize
#' @examples

runPhylo <- function(scCNA,
                     method = "nj",
                     metric = "euclidean",
                     n_threads = 1) {

  # checking distance matrix
  if (length(copykit::distMat(scCNA)) == 0) {
    message("No distance matrix detected in the scCNA object.")
    scCNA <-  runDistMat(scCNA, metric = metric)
  }

  if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
    stop("Number of samples in the distance matrix different from number of samples in the scCNA object. Perhaps you filtered your dataset? use copykit::runDistMat() to update it.")
  }

  # getting data
  seg_data <- t(segment_ratios(scCNA)) %>%
    as.data.frame()

  # ordering cells
  message("Creating neighbor-joining tree.")
  tree <- ape::nj(distMat(scCNA))

  tree <- ape::ladderize(tree)

  phylo(scCNA) <- tree

  message("Access slot with copykit::phylo(scCNA).")
  message("Done.")
  return(scCNA)

}

