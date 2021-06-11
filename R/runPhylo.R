#' Run phylogenetic analysis
#'
#' Performs phylogenetic analysis
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param method Phylogenetic method to be run, currently accepts "nj" (neighbor-joining) and "me" (minimum evolution). Defaults to "nj".
#' @param metric distance metric passed to construct the phylogeny (Defaults to "euclidean").
#' @param integer Whether the analysis is performed on the integer data. Used with integer_slot. Defaults to FALSE (using segment ratios data).
#' @param integer_slot Assay name where integer data is saved.
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
                     integer =  FALSE,
                     integer_slot,
                     n_threads = parallel::detectCores() / 4) {

  # cores check
  if (n_threads < 1) {
    n_threads <- 1
  }


  # getting data

  ## with ratios
  if (! integer) {

    message("Using ratio data...")
    seg_data <- segment_ratios(scCNA)
    seg_data[, ncol(seg_data)+1] <- 1
    seg_data[, ncol(seg_data)+1] <- 1
    seg_data <- t(seg_data) %>% as.data.frame()

  }  

  ## with integers
  if (integer) {

    if (is.null(integer_slot)){

      stop("Please specifiy integer slot.")

    } else if ( (integer_slot %in% names(SummarizedExperiment::assays(scCNA))) & 
                (integer_slot %!in% c("segment_ratios", "ratios", "bin_counts")) ) {

      message(sprintf("Using %s data...", integer_slot))
      seg_data <- SummarizedExperiment::assay(scCNA, integer_slot)
      seg_data[, ncol(seg_data)+1] <- 2
      seg_data[, ncol(seg_data)+1] <- 2
      seg_data <- t(seg_data) %>% as.data.frame()

    } else {

      stop("No integer data found in the integer slot! Please check the assays.")

    }

  }

  # calculating distance matrix
  message("Calculating distance matrix")
  distMat <- amap::Dist(seg_data,
                        method = metric,
                        nbproc = n_threads)

  # ordering cells
  if ( method %in% c("nj","me") ) {
    if ( method == "nj" ) {
       message("Creating neighbor-joining tree.")
       tree <- ape::nj(distMat)
    }

    if ( method == "me" ) {
       message("Creating minimum evolution tree.")
       tree <- ape::fastme.bal(distMat)
    }

  } else {
    stop("Currently only nj and me trees are supported.")
  }
  

  # root the tree
  tree <- ape::root.phylo(tree,
             outgroup = which(tree$tip.label == paste0("V",Ntip(tree))),
             resolve.root = T)
  tree <- ape::drop.tip(tree, tip = as.character(c(
              paste0("V",nrow(seg_data)), paste0("V",nrow(seg_data) - 1)
            )))

  tree <- ape::ladderize(tree)

  phylo(scCNA) <- tree

  message("Access slot with copykit::phylo(scCNA).")
  message("Done.")
  return(scCNA)

}

