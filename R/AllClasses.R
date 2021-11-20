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
#'  @importMethodsFrom SummarizedExperiment colData
#'  @importMethodsFrom SingleCellExperiment SingleCellExperiment
#'  @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SingleCellExperiment
#'  @importClassesFrom S4Vectors DataFrame SimpleList

#' A S4 class to store copy number assays with CopyKit.
#' Inherits from SingleCellExperiment
#'
#' @slot phylo Stores the single cell phylogenetic information with ape class
#' phylo
#' @slot consensusPhylo Stores the consensus phylogenetic information with
#' ape class phylo
#' @slot distMat Stores a distance matrix object used for graphs and heatmaps
#' @slot graph Stores an igraph object for network based clustering
#' @slot consensus stores a consensus data frame from \link{calcConsensus}
#' @export
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
.CopyKit <- setClass(
  "CopyKit",
  slots = representation(
    phylo = "phylo",
    consensusPhylo = "phylo",
    distMat = "dist",
    graph = "igraph",
    consensus = "data.frame"
  ),
  contains = "SingleCellExperiment"
)

#' @export
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
CopyKit <- function(consensus = data.frame(),
                    phylo = structure(list(), class = "phylo"),
                    consensusPhylo = structure(list(), class = "phylo"),
                    distMat = dist(matrix(0, 0, 0)),
                    graph = igraph::graph.empty(),
                    ...) {
  cna <-
    SingleCellExperiment::SingleCellExperiment(...)
  .CopyKit(
    cna,
    phylo = phylo,
    consensusPhylo = consensusPhylo,
    distMat = distMat,
    graph = graph,
    consensus = consensus
  )
}
