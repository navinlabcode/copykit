#' Find Clusters
#'
#' Search for clusters in the scCNA data.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param embedding String with the name of the reducedDim to pull data from.
#' @param ncomponents An integer with the number of components dimensions to
#' use from the embedding.
#' @param method A string with method used for clustering.
#' @param k_superclones A numeric k-nearest-neighbor value.
#' Used to find the superclones.
#' @param k_subclones A numeric k-nearest-neighbor value.
#' Used to find the subclones
#' @param seed A numeric passed on to pseudo-random dependent functions.
#'
#' @details \code{findClusters} uses the reduced dimensional embedding resulting
#'  from \code{\link{runUmap}} to perform clustering at two levels, hereby
#'  referred to as superclones, and subclones. When clustering for superclones
#'  findClusters creates a graph representation of the dataset reduced
#'  dimension embedding using a shared nearest neighbor algorithm
#'  (SNN) \code{\link[scran]{buildSNNGraph}}, from this graph the connected
#'  components are extracted and generally represent high-level structures
#'  that share large, lineage defining copy number events. At a more
#'  fine-grained resolution, CopyKit can also be used to detect subclones,
#'  i. e. groups of cells containing a unique copy number event per cluster,
#'  to do so the umap embedding is again used as the pre-processing step,
#'  this time to perform a density-based clustering with hdbscan
#'   \code{\link[dbscan]{hdbscan}}. Network clustering
#'  algorithms on top of the SNN graph such as the leiden algorithm
#'  \code{\link[leidenbase]{leiden_find_partition}}.
#'
#'  \itemize{
#'  \item{hdbscan}: hdbscan is an outlier aware clustering algorithm, since
#'  extensive filtering of the dataset can be applied before clustering with
#'  \code{\link{findOutliers}}, any cell classified as an outlier is inferred
#'  to the same cluster group as its closest, non-outlier, nearest-neighbor
#'   according to Euclidean distance.
#'  }
#'
#' @return Cluster information is added to
#' \code{\link[SummarizedExperiment]{colData}} in columns superclones or
#' subclones. Superclones are prefixed by 's' whereas subclones are prefixed
#' by 'c'.
#'
#' @seealso \code{\link{findSuggestedK}}.
#'
#' @references Laks, E., McPherson, A., Zahn, H., et al. (2019). Clonal
#'  Decomposition and DNA Replication States Defined by Scaled Single-Cell
#'  Genome Sequencing. Cell, 179(5), 1207–1221.e22.
#'  https://doi.org/10.1016/j.cell.2019.10.026
#'
#' Leland McInnes and John Healy and James Melville. UMAP: Uniform Manifold
#' Approximation and Projection for Dimension Reduction. arXiv:1802.03426
#'
#' Lun ATL, McCarthy DJ, Marioni JC (2016). “A step-by-step workflow for
#' low-level analysis of single-cell RNA-seq data with Bioconductor.”
#' F1000Res., 5, 2122. doi: 10.12688/f1000research.9501.2.
#'
#' @seealso \code{\link[dbscan]{hdbscan}} For hdbscan clustering.
#'
#' @export
#' @importFrom tidyr gather
#' @importFrom dplyr filter slice_min group_by right_join ungroup
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom SummarizedExperiment colData
#' @importFrom scran buildSNNGraph
#' @importFrom dbscan hdbscan
#' @importFrom tibble rownames_to_column
#' @importFrom forcats fct_reorder
#' @importFrom gtools mixedsort
#' @importFrom igraph cluster_leiden membership cluster_louvain
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- findClusters(copykit_obj)
findClusters <- function(scCNA,
    embedding = "umap",
    ncomponents = 2,
    method = c("hdbscan", "leiden", "louvain"),
    k_superclones = NULL,
    k_subclones = NULL,
    seed = 17) {
    method <- match.arg(method)

    # Obtain reduced dim embedding subset by the number of components.
    red_dim <- as.data.frame(reducedDim(scCNA, embedding)[,1:ncomponents])

    # checks
    if (is.null(k_subclones) &&
        is.null(S4Vectors::metadata(scCNA)$suggestedK)) {
        stop("k_subclones must have a numeric value.")
    }

    # if suggestedK is not null, use it as default
    if (is.null(k_subclones) &&
        !is.null(S4Vectors::metadata(scCNA)$suggestedK)) {
        message(
            "Using suggested k_subclones = ",
            S4Vectors::metadata(scCNA)$suggestedK
        )

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
            scran::buildSNNGraph(red_dim, k = k_superclones, transposed = TRUE)
        superclones <-
            as.factor(paste0("s", igraph::membership(igraph::components(g_major))))
        # storing info
        SummarizedExperiment::colData(scCNA)$superclones <- superclones
    }

    message("Finding clusters, using method: ", method)

    # subclones using leiden
    if (method == "leiden" || method == "louvain") {
        g_minor <-
            scran::buildSNNGraph(red_dim, k = k_subclones, transposed = TRUE)

        # saving g_minor graph
        copykit::graph(scCNA) <- g_minor

        # subclones
        if (method == "leiden") {
            leid_obj <- igraph::cluster_leiden(g_minor,
                resolution_parameter = 0.2,
                n_iterations = 100
            )
        }

        if (method == "louvain") {
            leid_obj <- igraph::cluster_leiden(g_minor)
        }

        leid_obj_com <- igraph::membership(leid_obj)

        subclones <- as.factor(paste0("c", leid_obj$membership))

        n_clones <- length(unique(subclones))
        message("Found ", n_clones, " subclones.")
    }

    # subclones using hdbscan
    if (method == "hdbscan") {
        withr::with_seed(seed,
                         hdb <- dbscan::hdbscan(red_dim,
                                                minPts = k_subclones))

        hdb_clusters <- as.character(hdb$cluster)

        hdb_df <- data.frame(
            cell = rownames(red_dim),
            hdb = hdb_clusters
        )

        subclones <- as.factor(paste0("c", hdb_df$hdb))

        n_clones <- length(unique(subclones))
        n_outliers <- length(subclones[subclones == "c0"])

        if (n_outliers > 0) {
            message("Found ", n_outliers, " outliers cells (group 'c0')")
            message("Found ", n_clones - 1, " subclones.")
        } else {
            message("Found ", n_clones, " subclones.")
        }
    }

    # storing subclones info
    SummarizedExperiment::colData(scCNA)$subclones <-
        forcats::fct_relevel(
            droplevels(subclones),
            gtools::mixedsort(unique(as.character(subclones)))
        )

    message("Done.")

    return(scCNA)
}
