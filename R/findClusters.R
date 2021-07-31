#' Find Clusters
#'
#' Search for clusters in the scCNA data.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param embedding String with the name of the reducedDim to pull data from.
#' @param method A string with method used for clustering.
#' @param k_superclones A numeric scalar k-nearest-neighbor value.
#' Used to find the superclones.
#' @param k_subclones A numeric scalar k-nearest-neighbor value.
#' Used to find the subclones
#' @param seed A numeric scalar seed passed on to pseudo-random dependent functions.
#'
#' @details \code{findClusters} uses the reduced dimensional embedding resulting
#'  from \code{\link{runUmap}} to perform clustering at two levels, hereby referred
#'  to as superclones, and subclones. When clustering for superclones findClusters
#'  creates a graph representation of the dataset reduced dimension embedding
#'  using a shared nearest neighbor algorithm (SNN) \code{\link[scran]{buildSNNGraph}},
#'  from this graph the connected components are extracted and generally
#'  represent high-level structures that share large, lineage defining copy
#'  number events. At a more fine-grained resolution, CopyKit can also be
#'  used to detect subclones, i. e. groups of cells containing a unique
#'  copy number event per cluster, to do so the umap embedding is again
#'  used as the pre-processing step, this time to perform a density-based
#'  clustering with hdbscan \code{\link[dbscan]{hdbscan}}. Network clustering
#'  algorithms on top of the SNN graph such as the leiden algorithm
#'  \code{\link[leidenbase]{leiden_find_partition}}.
#'
#'  \itemize{
#'  \item{hdbscan}: hdbscan is an outlier aware clustering algorithm, since
#'  extensive filtering of the dataset can be applied before clustering with
#'  \code{\link{filterCells}}, any cell classified as an outlier is inferred
#'  to the same cluster group as its closest, non-outlier, nearest-neighbor
#'   according to Euclidean distance.
#'  }
#'
#' @return Cluster information is added to \code{\link[SummarizedExperiment]{colData}}
#' in columns superclones or subclones. Superclones are prefixed by 's' whereas subclones
#' are prefixed by 'c'
#'
#' @seealso \code{\link{findSuggestedK}} to obtain suggestions of k_subclones values.
#'
#' @references Laks, E., McPherson, A., Zahn, H., et al. (2019). Clonal Decomposition
#' and DNA Replication States Defined by Scaled Single-Cell Genome Sequencing.
#' Cell, 179(5), 1207–1221.e22. https://doi.org/10.1016/j.cell.2019.10.026
#'
#' Leland McInnes and John Healy and James Melville. UMAP: Uniform Manifold
#' Approximation and Projection for Dimension Reduction. arXiv:1802.03426
#'
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for low-level
#' analysis of single-cell RNA-seq data with Bioconductor.”
#' F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#'
#' @seealso \code{\link[dbscan]{hdbscan}} For hdbscan clustering.
#'
#' @export
#' @import leidenbase
#' @importFrom tidyr gather
#' @importFrom dplyr filter slice_min group_by right_join ungroup
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @importFrom scran buildSNNGraph
#' @importFrom dbscan hdbscan
#' @importFrom tibble rownames_to_column
#' @importFrom forcats fct_reorder
#' @importFrom gtools mixedsort
#'
#' @examples

findClusters <- function(scCNA,
                         embedding = "umap",
                         method = c("hdbscan", "leiden"),
                         k_superclones = NULL,
                         k_subclones = NULL,
                         seed = 17) {

  method <- match.arg(method)

  # obtaining data from reducedDim slot
  if (!is.null(SingleCellExperiment::reducedDim(scCNA, embedding))) {

    umap_df <-
      SingleCellExperiment::reducedDim(scCNA, embedding) %>%
      as.data.frame()

  } else
    stop("Reduced dimensions slot is NULL. Use runUmap() to create it.")

  # checks
  if (is.null(k_subclones) && is.null(S4Vectors::metadata(scCNA)$suggestedK)) {
    stop("k_subclones must have a numeric value.")
  }

  # if suggestedK is not null, use it as default
  if (is.null(k_subclones) && !is.null(S4Vectors::metadata(scCNA)$suggestedK)) {
    message(paste("Using suggested k_subclones =",
                  S4Vectors::metadata(scCNA)$suggestedK))

    k_subclones <- S4Vectors::metadata(scCNA)$suggestedK
  }

  if (!is.numeric(k_subclones)) {
    stop("k_subclones must be a numeric values.")
  }

  # superclones clustering
  if (!is.null(k_superclones)) {

    # type check
    if (!is.numeric(k_superclones)) {
      stop("k_superclones must have a numeric value.")
    }

    g_major <-
      scran::buildSNNGraph(umap_df, k = k_superclones, transposed = T)
    superclones <-
      as.factor(paste0("s", igraph::membership(igraph::components(g_major))))
    #storing info
    SummarizedExperiment::colData(scCNA)$superclones <- superclones

  }

  message(paste("Finding clusters, using method:", method))

  # subclones using leiden
  if (method == "leiden") {

    g_minor  <-
      scran::buildSNNGraph(umap_df, k = k_subclones, transposed = T)
    g_adj <- igraph::as_adjacency_matrix(g_minor)

    # saving g_minor graph
    copykit::graph(scCNA) <- g_minor

    #minor
    leid_obj <- try(leidenbase::leiden_find_partition(
      g_minor,
      partition_type = 'RBConfigurationVertexPartition',
      resolution_parameter = 1,
      seed = seed
    ))
    if (inherits(leid_obj, "try-error")) {
      stop('Running leiden failed.')
    } else {
      subclones <- as.factor(paste0('c', leid_obj$membership))
    }

    n_clones <- length(unique(subclones))
    message(paste("Found", n_clones, "subclones."))

  }

  # subclones using hdbscan
  if (method == "hdbscan") {
    set.seed(seed)
    hdb <- dbscan::hdbscan(umap_df,
                           minPts = k_subclones)
    hdb_clusters <- as.character(hdb$cluster)

    hdb_df <- data.frame(cell = rownames(umap_df),
                         hdb = hdb_clusters)

    subclones <- as.factor(paste0('c', hdb_df$hdb))

    n_clones <- length(unique(subclones))
    n_outliers <- length(subclones[subclones == 'c0'])

    message(paste("Found", n_clones, "subclones."))
    message(paste(n_outliers,
                  "cells were classified as outliers. Check subclone group 'c0'."))


  }

  # storing subclones info
  SummarizedExperiment::colData(scCNA)$subclones <-
    forcats::fct_relevel(droplevels(subclones),
                         gtools::mixedsort(unique(as.character(subclones))))

  message("Done.")

  return(scCNA)

}
