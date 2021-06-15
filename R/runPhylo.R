#' Run phylogenetic analysis
#'
#' Performs phylogenetic analysis
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param method Phylogenetic method to be run, currently accepts "nj" (neighbor-joining) and "me" (minimum evolution). Defaults to "nj".
#' @param metric distance metric passed to construct the phylogeny (Defaults to "euclidean").
#' @param integer Whether the analysis is performed on the integer data. Defaults to FALSE (using segment ratios data).
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `amap::Dist`
#'
#' @return A rooted phylogenetic tree object in the slot \code{phylo} from scCNA object. Access phylo slot with: \code{copykit::phylo(scCNA)}
#' @export
#'
#' @importFrom ape nj fastme.bal ladderize
#' @examples

runPhylo <- function(scCNA,
                     method = "nj",
                     metric = "euclidean",
                     integer =  FALSE,
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

    if ("integer" %in% names(SummarizedExperiment::assays(scCNA))) {

      message("Using integer data...")
      seg_data <- SummarizedExperiment::assay(scCNA, "integer")
      seg_data[, ncol(seg_data)+1] <- 2
      seg_data[, ncol(seg_data)+1] <- 2
      seg_data <- t(seg_data) %>% as.data.frame()

    } else {

      stop("Integer data not found! Please make sure the integer data is in SummarizedExperiment::assay(scCNA, \"integer\").")

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

