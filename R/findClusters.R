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
#' @param method Which method should be used for clustering, options are "hdbscan" or "leiden". Defaults to "hdbscan".
#' @param k_superclones k-nearest-neighbor value that will be used to find the major clusters.(Defaults to 35).
#' @param k_subclones k-nearest-neighbor value that will be used to find the minor clusters (Defaults to 21).
#' @param seed Seed passed on to leiden algorithm (Defaults to 17).
#'
#' @return Metadata cluster information that can be found in \code{SummarizedExperiment::colData(scCNA)$superclones} for the major clusters and \code{SummarizedExperiment::colData(scCNA)$subclones} for the minor clusters. Major clusters are named with capital letters whereas minor clusters are named with a numbers but as a character vector.
#'
#' @export
#' @import leidenbase
#' @importFrom tidyr gather
#' @importFrom dplyr filter slice_min group_by right_join
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @importFrom scran buildSNNGraph
#' @importFrom dbscan hdbscan
#' @importFrom tibble rownames_to_column
#'
#' @examples

findClusters <- function(scCNA,
                         method = "hdbscan",
                         k_superclones = NULL,
                         k_subclones = NULL,
                         seed = 17) {

  # obtaining data from reducedDim slot
  if (!is.null(SingleCellExperiment::reducedDim(scCNA))) {
    umap_df <-
      SingleCellExperiment::reducedDim(scCNA, 'umap') %>%
      as.data.frame()

  } else
    stop("Reduced dimensions slot is NULL. Use runUmap() to create it.")

  # settting k defaults
  if (is.null(k_superclones)) {
    k_superclones <- 0.05*nrow(umap_df)
  }

  if (is.null(k_subclones)) {
    k_subclones <- 0.02*nrow(umap_df)
  }

  # checks
  if (!is.numeric(k_subclones) || !is.numeric(k_superclones)) {
    stop("k_subclones and k_superclones must be numeric values.")
  }

  g_major <-
    scran::buildSNNGraph(umap_df, k = k_superclones, transposed = T)
  g_minor  <-
    scran::buildSNNGraph(umap_df, k = k_subclones, transposed = T)
  g_adj <- igraph::as_adjacency_matrix(g_minor)

  # saving g_minor graph
  copykit::graph(scCNA) <- g_minor

  # finding clusters
  #  major
  message(paste("Finding clusters, using method:", method))
  superclones <- as.factor(paste0("s", igraph::membership(igraph::components(g_major))))

  # using leiden
  if (method == "leiden") {
    #minor
    leid_obj <- try(leidenbase::leiden_find_partition(
      g_minor,
      partition_type = 'RBConfigurationVertexPartition',
      resolution_parameter = 1,
      seed=seed))
    if (inherits(leid_obj, "try-error")) {
      stop('Running leiden failed.')
    } else {
      subclones <- as.factor(paste0('c', leid_obj$membership))
    }
  }

  # using hdbscan
  if (method == "hdbscan") {

    set.seed(seed)
    hdb <- dbscan::hdbscan(umap_df,
                           minPts = k_subclones)
    hdb_clusters <- as.character(hdb$cluster)

    # hdbscan is an outlier aware clustering algorithm. Copykit assumes that filterCells already took care of removing bad cells
    # this is why any outlier is added to the ones classified as outliers to the closest cluster possible according to euclidean distance

    dist_umap <- dist(umap_df) %>%
      as.matrix() %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell2") %>%
      tidyr::gather(key = "cell1",
                    value = "dist",
                    -cell2) %>%
      dplyr::filter(cell1 != cell2)

    hdb_df <- data.frame(cell = rownames(umap_df),
                         hdb = hdb_clusters)

    dist_min <- dist_umap %>%
      dplyr::right_join(hdb_df, by = c("cell2" = "cell")) %>%
      dplyr::filter(hdb != "0") %>%
      dplyr::group_by(cell1) %>%
      dplyr::slice_min(dist) %>%
      ungroup()

    for (i in 1:nrow(umap_df)) {

      if(hdb_df$hdb[i] == "0") {
        cellname <- rownames(umap_df)[i]
        closest_cell <- dplyr::filter(dist_min, cell1 == rownames(umap_df)[i])$cell2
        closest_cell_cluster <- dplyr::filter(hdb_df, cell == closest_cell)$hdb
        hdb_df$hdb[i] <- closest_cell_cluster
      }

      subclones <- as.factor(paste0('c',hdb_df$hdb))

    }

  }

  # storing info
  SummarizedExperiment::colData(scCNA)$superclones <- superclones
  SummarizedExperiment::colData(scCNA)$subclones <- droplevels(subclones)

  message("Done.")

  return(scCNA)

}


