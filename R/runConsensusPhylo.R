#' runConsensusPhylo
#'
#'  Runs a minimal evolution tree algorithm for the consensus data frame
#'
#' @param consensus_df The consensus data frame from the segmented dataset
#' @param clusters The results from clustering
#' @param ploidy_VAL The inferred ploidy value. Will be used to scale the data
#' @param rotate_nodes Nodes to be rotated
#' @param plot Boolean. Prints the plot if true
#'
#' @importFrom ape fastme.bal Ntip root.phylo drop.tip
#'
#' @return
#' @export
#'
#' @examples
runConsensusPhylo <- function(scCNA) {

  if (nrow(consensus(scCNA)) == 0) {
    stop("Consensus slot is empty. run calcConsensus().")
  }

  consensus_df <- as.data.frame(t(consensus(scCNA)))

  #adding a neutral state, will use as root
  consensus_df[nrow(consensus_df) + 1, ] <- 1
  consensus_df[nrow(consensus_df) + 1, ] <- 1

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Nov 20 12:24:27 2020
  # tree ME
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Nov 20 12:24:35 2020

  tree <- ape::fastme.bal(dist(consensus_df, method = "manhattan"))

  tree <-
    ape::root.phylo(tree,
               outgroup = which(tree$tip.label == ape::Ntip(tree)),
               resolve.root = T)

  tree <-
    ape::drop.tip(tree, tip = as.character(c(
      nrow(consensus_df), nrow(consensus_df) - 1
    )))

  tree <- ladderize(tree)

  consensusPhylo(scCNA) <- tree

  return(scCNA)

}
