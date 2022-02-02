### all clases for the copykit package

###################################################################
# ape has no formal definition of the phylo-class
# that can be used for @importClassesFrom
# the solution adopted here is the same as the one adopted phyloseq package
# which creates S3 and S4 placeholders for the phylo class

#' @keywords internal
phylo <- structure(list(), class = "phylo")

#' Placeholder for the igraph class
#' @exportClass igraph
#' @name igraph-class
#' @rdname CopyKit-class
#' @keywords internal
setOldClass("igraph")

###################################################################
#' The CopyKit class
#'
#' S4 Class that extends the Bioconductor SingleCellExperiment class to hold
#' single cell copy number datasets.
#'
#' @slot phylo Stores the single cell phylogenetic information with ape class
#' phylo.
#' @slot consensusPhylo Stores the consensus phylogenetic information with
#' ape class phylo.
#' @slot distMat Stores a distance matrix object used for graphs and heatmaps.
#' @slot graph Stores an igraph object for network based clustering.
#' @slot consensus stores a consensus data frame from
#' \code{\link{calcConsensus}.}
#' @return A CopyKit class object.
#' @references The Bioconductor SingleCellExperiment Class
#' DOI: 10.18129/B9.bioc.SingleCellExperiment
#' @import methods
#' @name CopyKit-class
#' @rdname CopyKit-class
#' @importMethodsFrom SummarizedExperiment colData
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
#' @exportClass CopyKit
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
#' @rdname CopyKit-class
#' @param consensus A data frame with the consensus information.
#' @param phylo A phylo object with a phylogenetic tree.
#' @param consensusPhylo A phylo object with a phylogenetic consensus tree.
#' @param graph A graph object with a graph made from the umap data.
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
