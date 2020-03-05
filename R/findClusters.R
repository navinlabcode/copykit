#' Find Clusters
#'
#' Search for clusters in the scCNA data by using a graph based approach. \code{findClusters()} builds an SNN graph of the k-nearest neighbors and attempts to find two different configuration of clusters. Major and minor subpopulations.
#' Major clusters are found by looking at the graph connected components, whereas the minor clusters use the Leiden algorithm to detect connected communities within the major clusters.
#' \code{findClusters()} generates the graph by using the UMAP embedding that can be obtained after running \code{runUmap()}.
#'
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param k_major k-nearest-neighbor value that will be used to find the major clusters. Must be higher than k_minor (Defaults to 35).
#' @param k_minor k-nearest-neighbor value that will be used to find the minor clusters (Defaults to 21).
#' @param seed Seed passed on to leiden algorithm (Defaults to 17).
#'
#' @return Metadata cluster information that can be found in \code{SummarizedExperiment::colData(scCNA)$major_clusters} for the major clusters and \code{SummarizedExperiment::colData(scCNA)$minor_clusters} for the minor clusters. Major clusters are named with capital letters whereas minor clusters are named with a numbers but as a character vector.
#'
#' @import leidenbase
#'
#' @export
#'
#' @examples

findClusters <- function(scCNA,
                         k_major = 35,
                         k_minor = 21,
                         seed = 17) {
  # checks
  if (!is.numeric(k_minor) || !is.numeric(k_major)) {
    stop("k_minor and k_major must be numeric values.")
  }

  # obtaining data from reducedDim slot
  if (!is.null(SingleCellExperiment::reducedDim(scCNA))) {
    umap_df <-
      SingleCellExperiment::reducedDim(scCNA, 'umap') %>%
      as.data.frame()

  } else
    stop("Reduced dimensions slot is NULL. Use runUmap() to create it.")

  # building graph
  message("Building SNN graph.")
  if (k_major < k_minor) {
    stop("k_major argument needs to be equal or higher than k_minor argument.")
  }

  g_major <-
    scran::buildSNNGraph(umap_df, k = k_major, transposed = T)
  g_minor  <-
    scran::buildSNNGraph(umap_df, k = k_minor, transposed = T)
  g_adj <- igraph::as_adjacency_matrix(g_minor)

  # saving g_minor graph
  copykit::graph(scCNA) <- g_minor

  # finding clusters
  #  major
  message("Finding clusters.")
  g_clusters <- igraph::membership(igraph::components(g_major))
  g_clusters <-
    sapply(strsplit(paste(g_clusters), ''), function(y)
      paste(LETTERS[as.numeric(y)], collapse = ''))

  #minor
  # leid <- try(leiden::leiden(g_adj, seed = seed))
  leid_obj <- try(leiden_find_partition(
    g_minor,
    partition_type = 'RBConfigurationVertexPartition',
    resolution_parameter = 1))
  if (inherits(leid_obj, "try-error")) {
    leid <- g_clusters
    warning('Running leiden fails. Copy major clusters to minor clusters.')
  } else {
    leid <- leid_obj$membership
  }

  # storing info
  SummarizedExperiment::colData(scCNA)$major_clusters <- g_clusters
  SummarizedExperiment::colData(scCNA)$minor_clusters <- leid

  message("Done.")

  return(scCNA)

}
