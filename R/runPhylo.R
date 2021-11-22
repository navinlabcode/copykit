#' Run phylogenetic analysis
#'
#' Performs phylogenetic analysis
#'
#' @author Darlan Conterno Minussi
#' @author Junke Wang
#'
#' @param scCNA scCNA object.
#' @param method Phylogenetic method to be run, currently accepts "nj" (neighbor-joining) and "me" (minimum evolution). Defaults to "nj".
#' @param metric distance metric passed to construct the phylogeny (Defaults to "euclidean").
#' @param assay String with the name of the assay to pull data from to run phylogenetic analysis. Note that only assay named "integer" will be treated as integer.
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `amap::Dist`
#'
#' @return A rooted phylogenetic tree object in the slot \code{phylo} from scCNA object. Access phylo slot with: \code{copykit::phylo(scCNA)}
#' @export
#'
#' @importFrom ape nj fastme.bal ladderize
#' @examples
#' copykit_obj <- copykit_example()
#' copykit_obj <- findNormalCells(copykit_obj)
#' copykit_obj <- copykit_obj[,colData(copykit_obj)$is_normal == "FALSE"]
#' copykit_obj <- filterCells(copykit_obj)
#' copykit_obj <- copykit_obj[,colData(copykit_obj)$filtered == "kept"]
#' copykit_obj <- runPhylo(copykit_obj)
#'

runPhylo <- function(scCNA,
                     method = "nj",
                     metric = "euclidean",
                     assay = "segment_ratios",
                     n_threads = parallel::detectCores() / 4) {
  # cores check
  if (n_threads < 1) {
    n_threads <- 1
  }


  # getting data
  if (!assay %in% names(SummarizedExperiment::assays(scCNA))) {
    stop("No data found in the assay! Please check the assay name.")
  }

  seg_data <- SummarizedExperiment::assay(scCNA, assay)


  if (assay == "integer") {
    ## with integers
    message("Using integer data...")
    seg_data[, ncol(seg_data) + 1] <- 2
    seg_data[, ncol(seg_data) + 1] <- 2
    seg_data <- t(seg_data) %>% as.data.frame()

  } else {
    # with ratios
    message("Using ratio data...")
    seg_data[, ncol(seg_data) + 1] <- 1
    seg_data[, ncol(seg_data) + 1] <- 1
    seg_data <- t(seg_data) %>% as.data.frame()

  }


  # calculating distance matrix
  message("Calculating distance matrix")
  distMat <- amap::Dist(seg_data,
                        method = metric,
                        nbproc = n_threads)

  # ordering cells
  if (method %in% c("nj", "me")) {
    if (method == "nj") {
      message("Creating neighbor-joining tree.")
      tree <- ape::nj(distMat)
    }

    if (method == "me") {
      message("Creating minimum evolution tree.")
      tree <- ape::fastme.bal(distMat)
    }

  } else {
    stop("Currently only nj and me trees are supported.")
  }


  # root the tree
  tree <- ape::root.phylo(tree,
                          outgroup = which(tree$tip.label == paste0("V", Ntip(tree))),
                          resolve.root = TRUE)
  tree <- ape::drop.tip(tree, tip = as.character(c(
    paste0("V", nrow(seg_data)), paste0("V", nrow(seg_data) - 1)
  )))

  tree <- ape::ladderize(tree)

  phylo(scCNA) <- tree

  message("Access slot with copykit::phylo(scCNA).")
  message("Done.")
  return(scCNA)

}
