### all clases for the copykit package

###################################################################
# ape has no formal definition of the phylo-class
# that can be used for @importClassesFrom
# the solution adopted here is the same as the one adopted phyloseq package
# which creates S3 and S4 placeholders for the phylo class

#' @keywords internal
phylo <- structure(list(), class = "phylo")

#' @exportClass igraph
setOldClass("igraph")

###################################################################
# defining scCNA class
#'
#'  @export
#'  @importMethodsFrom SingleCellExperiment SingleCellExperiment
#'  @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SingleCellExperiment
#'  @importClassesFrom S4Vectors DataFrame SimpleList


#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.scCNA <- setClass("scCNA",
                   slots = representation(
                     phylo = "phylo",
                     consensusPhylo = "phylo",
                     distMat = "dist",
                     graph = "igraph",
                     consensus = "data.frame"
                   ),
                   contains = "SingleCellExperiment")

#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
scCNA <- function(segment_ratios,
                  ratios,
                  bin_counts,
                  consensus = data.frame(),
                  phylo = structure(list(), class = "phylo"),
                  consensusPhylo = structure(list(), class = "phylo"),
                  distMat = dist(matrix(0,0,0)),
                  graph = igraph::graph.empty(),
                  ...) {
  cna <-
    SingleCellExperiment::SingleCellExperiment(list(segment_ratios = segment_ratios,
                                                    ratios = ratios,
                                                    bin_counts = bin_counts),
                                               ...)
  .scCNA(cna,
         phylo = phylo,
         consensusPhylo = consensusPhylo,
         distMat = distMat,
         graph = graph,
         consensus = consensus)
}




