#' runConsensusPhylo
#'
#'  Runs a minimal evolution tree algorithm for the consensus data frame
#'
#' @param scCNA The scCNA object.
#' @param root A string indicating how to root the consensus tree.
#'
#' @importFrom ape fastme.bal Ntip root.phylo drop.tip
#'
#' @return
#' @export
#'
#' @examples
runConsensusPhylo <- function(scCNA,
                              root = c('mrca', 'neutral')) {

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
