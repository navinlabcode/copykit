#' findSuggestedK
#'
#' Performs a grid search over a range of k values to assess cluster stability.
#'
#' @param scCNA  scCNA object.
#' @param embedding String with the name of the reducedDim embedding to pull data from.
#' @param k_range A numeric range of values to be tested.
#' @param method A string with the method of clustering to be tested.
#' @param metric A string with the function to summarize the jaccard similarity
#' value from all clusters.
#' @param seed A numerical scalar with a seed value to be passed on to
#' \code{\link[uwot]{umap}}.
#' @param B A numeric with the number of bootstrapping iterations passed on to
#' \code{\link[fpc]{clusterboot}}. Higher values yield better results at a cost
#' of performance
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details Performs a grid-search over a range of k values and returns the value
#' that maximizes the jaccard similarity. Importantly, while this approach does
#' not guarantee optimal clustering, it provides a guide that maximizes cluster
#' stability.
#'
#' The default tested range is from 7 to the square root of the number of cells
#' in the scCNA object. If sqrt(n_cells) is smaller than 7 a range of 5 to 15
#' is tested.
#'
#' @return Adds a table with the mean jaccard coefficient of clusters for each
#' tested k and the suggested k value to be used for clustering to
#' \code{\link[SummarizedExperiment]{metadata}}
#'
#' @seealso \code{\link[fpc]{clusterboot}}
#' @seealso \code{\link{plotSuggestedK}}
#'
#' @references Hennig, C. (2007) Cluster-wise assessment of cluster stability.
#' Computational Statistics and Data Analysis, 52, 258-271.
#'
#' Hennig, C. (2008) Dissolution point and isolation robustness: robustness
#' criteria for general cluster analysis methods.
#' Journal of Multivariate Analysis 99, 1154-1176.
#'
#' @export
#'
#' @importFrom fpc clusterboot
#' @importFrom tibble enframe
#' @importFrom dplyr group_by summarise
#' @importFrom dbscan hdbscan
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment reducedDim
#'
#' @examples
findSuggestedK <- function(scCNA,
                           embedding = 'umap',
                           k_range = NULL,
                           method = c("hdbscan", "leiden"),
                           metric = c("mean", "median"),
                           seed = 17,
                           B = 200,
                           BPPARAM = bpparam())
{

  metric <- match.arg(metric)
  method <- match.arg(method)

  # defining k_range
  if (is.null(k_range)) {

    n_cells <- ncol(segment_ratios(scCNA))
    max_range <- round(sqrt(n_cells))
    min_range <- 5

    if (min_range > max_range) {
      k_range <- 5:20
    } else {
      k_range <- min_range:max_range
    }

  }

  # obtaining data from reducedDim slot
  if (is.null(SingleCellExperiment::reducedDim(scCNA, embedding))) {
    stop("Reduced dimensions slot is NULL. Use runUmap().")
  }

  message(cat("Calculating jaccard similarity for k range:", k_range))

  if (method == 'leiden') {

    # setting seed to a different variable to avoid recursive call
    seed_val <- seed

    jaccard <- BiocParallel::bplapply(k_range, function(i) {
      df_clusterboot <-
        .quiet(
          fpc::clusterboot(
            SingleCellExperiment::reducedDim(scCNA, "umap"),
            B = B,
            clustermethod = leidenCBI,
            seed_leid = seed_val,
            k = i
          )
        )

      n_subclones <- 1:df_clusterboot$nc
      n_cells_per_subclone <- as.numeric(table(df_clusterboot$partition))

      df_result <- data.frame(k = i,
                              subclones = paste0('c', n_subclones),
                              n_cells = n_cells_per_subclone,
                              bootmean = df_clusterboot$bootmean)

    }, BPPARAM = BPPARAM)

  }

  if (method == 'hdbscan') {

    jaccard <- BiocParallel::bplapply(k_range, function(i) {
      df_clusterboot <-
        .quiet(
          fpc::clusterboot(
            SingleCellExperiment::reducedDim(scCNA, "umap"),
            B = B,
            clustermethod = hdbscanCBI,
            seed = seed,
            minPts = i,
            noisemethod=TRUE
          )
        )

      # the outlier partition is the last element see ?fpc::clusterboot
      # here the oulier partition will be named as 'c0'
      n_subclones <- 1:df_clusterboot$nc
      names(n_subclones) <- paste0('c', n_subclones)
      n_cells_per_subclone <- as.numeric(table(df_clusterboot$partition))
      names(n_subclones)[length(n_subclones)] <- 'c0'

      df_result <- data.frame(k = i,
                              subclones = names(n_subclones),
                              n_cells = n_cells_per_subclone,
                              bootmean = df_clusterboot$bootmean)

    }, BPPARAM = BPPARAM)

  }

  # Obtaining values from the list and binding to a dataframe
  names(jaccard) <- k_range
  bootmean_df <- do.call(rbind, jaccard)
  jc_df <- bootmean_df %>%
    dplyr::group_by(k) %>%
    dplyr::summarise(mean = mean(bootmean),
              median = median(bootmean))


  # Using selected summarising function
  if (metric == 'mean') {

    jc_df_opt <- jc_df %>%
      dplyr::filter(mean == max(jc_df$mean))

    selected_k_jac_value <- jc_df_opt$mean

  }

  if (metric == 'median') {

    jc_df_opt <- jc_df %>%
      filter(median == max(jc_df$median))

    selected_k_jac_value <- jc_df_opt$median

  }

  # If the best jaccard similarity score is found in more than one k value
  # it prioritizes the maximum value. If the score remains tied (possible
  # when score == 1) prioritizes the smaller k
  if (nrow(jc_df_opt) > 1) {

    jc_df_opt <- jc_df_opt %>%
      filter(k == min(k))

    selected_k_jac_value <- jc_df_opt$median

  }

  message(paste(
    "Suggested k =",
    jc_df_opt$k,
    "with", metric, "jaccard similarity of:",
    round(selected_k_jac_value, 3)
  ))

  # adding values to metadata
  S4Vectors::metadata(scCNA)$suggestedK_df <- bootmean_df
  S4Vectors::metadata(scCNA)$suggestedK <- jc_df_opt$k

  return(scCNA)

}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# custom hdbscan function to pass to fpc::clusterboot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @internal
#' @export
hdbscanCBI <-
  function(data, minPts, diss = inherits(data, "dist"), ...) {

    if (diss)
      c1 <- dbscan::hdbscan(data, minPts, method = "dist", ...)
    else
      c1 <- dbscan::hdbscan(data, minPts, ...)

    partition <- c1$cluster
    cl <- list()
    nccl <- max(partition)
    partition[partition == 0] <- nccl + 1
    nc <- max(partition)
    for (i in 1:nc) cl[[i]] <- partition == i
    out <- list(result = c1, nc = nc, nccl = nccl, clusterlist = cl,
                partition = partition, clustermethod = "dbscan")

    out
  }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# custom leiden function to pass to fpc::clusterboot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' @internal
#' @export
leidenCBI <- function(data,k,seed_leid,diss=inherits(data,"dist"),...){

  g_minor  <-
    scran::buildSNNGraph(data, k = k, transposed = T)

  if (diss)
    c1 <- leidenbase::leiden_find_partition(
      g_minor,
      partition_type = 'RBConfigurationVertexPartition',
      resolution_parameter = 1,
      seed = seed_leid
    )
  else
    c1 <- leidenbase::leiden_find_partition(
      g_minor,
      partition_type = 'RBConfigurationVertexPartition',
      resolution_parameter = 1,
      seed = seed_leid
    )

  partition <- c1$membership

  cl <- list()
  nc <- max(partition)

  for (i in 1:nc)
    cl[[i]] <- partition==i
  out <- list(result=c1,nc=nc,clusterlist=cl,partition=partition,
              clustermethod="leiden")
  out
}

