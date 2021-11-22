#' runConsensusPhylo
#'
#'  Runs a minimal evolution tree algorithm for the consensus data frame
#'
#' @param scCNA The scCNA object.
#' @param root A string indicating how to root the consensus tree.
#' @param root_user A numeric with the vector to be used as root of the tree if
#' \code{root} is set to 'user'. Must have the same length as the number of bins
#' of the genome scaffold.
#'
#' @importFrom ape fastme.bal Ntip root.phylo drop.tip
#'
#' @return A phylo object with a consensus tree stored in the consensusPhylo slot
#' of the CopyKit object.
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example()
#' copykit_obj <- findNormalCells(copykit_obj)
#' copykit_obj <- copykit_obj[,colData(copykit_obj)$is_normal == "FALSE"]
#' copykit_obj <- filterCells(copykit_obj)
#' copykit_obj <- copykit_obj[,colData(copykit_obj)$filtered == "kept"]
#' copykit_obj <- runUmap(copykit_obj)
#' copykit_obj <- findSuggestedK(copykit_obj)
#' copykit_obj <- findClusters(copykit_obj)
#' copykit_obj <- calcConsensus(copykit_obj)
#' copykit_obj <- runConsensusPhylo(copykit_obj)
#'
runConsensusPhylo <- function(scCNA,
                              root = c('mrca', 'neutral', 'user'),
                              root_user = NULL) {

  root <- match.arg(root)

  if (nrow(consensus(scCNA)) == 0) {
    stop("Consensus slot is empty. run calcConsensus().")
  }

  consensus_df <- as.data.frame(t(consensus(scCNA)))

  if (root == 'neutral') {

    #adding a neutral state, will use as root
    consensus_df[nrow(consensus_df) + 1, ] <- 1
    consensus_df[nrow(consensus_df) + 1, ] <- 1

  }

  if (root == 'mrca') {
    #obtain number closest to the ground state for each
    anc_profile <- apply(consensus_df,
                         2,
                         function(x) x[which.min(abs(x-1))] )

    consensus_df[nrow(consensus_df)+1,] <- anc_profile
    consensus_df[nrow(consensus_df)+1,] <- anc_profile

  }

  if (root == 'user') {

    if (length(root_user) != ncol(consensus_df)) {
      stop("Length of root_user argument must be the same as nrow(scCNA).")
    }

    anc_profile <- root_user

    consensus_df[nrow(consensus_df)+1,] <- anc_profile
    consensus_df[nrow(consensus_df)+1,] <- anc_profile

  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Nov 20 12:24:27 2020
  # tree ME
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Nov 20 12:24:35 2020

  tree <- ape::fastme.bal(dist(consensus_df, method = "manhattan"))

  tree <-
    ape::root.phylo(tree,
               outgroup = which(tree$tip.label == ape::Ntip(tree)),
               resolve.root = TRUE)

  tree <-
    ape::drop.tip(tree, tip = as.character(c(
      nrow(consensus_df), nrow(consensus_df) - 1
    )))

  tree <- ladderize(tree)

  consensusPhylo(scCNA) <- tree

  return(scCNA)

}
