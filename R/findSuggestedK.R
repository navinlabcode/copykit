#' findSuggestedK
#'
#' Performs a grid search over a range of k values to assess cluster stability.
#'
#' @param scCNA  scCNA object.
#' @param embedding String with the name of the reducedDim embedding to pull data from.
#' @param k_range A numeric range of values to be tested.
#' @param method A string with the method of clustering to be tested.
#' @param fun A string with the function to summarize the jaccard similarity
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
                           method = "hdbscan",
                           fun = "median",
                           seed = 17,
                           B = 200,
                           BPPARAM = bpparam())
{

  # defining k_range
  if (is.null(k_range)) {

    n_cells <- ncol(segment_ratios(scCNA))
    max_range <- round(sqrt(n_cells))
    min_range <- 7

    if (min_range > max_range) {
      k_range <- 5:15
    } else {
      k_range <- min_range:max_range
    }

  }

  # obtaining data from reducedDim slot
  if (is.null(SingleCellExperiment::reducedDim(scCNA, embedding))) {
    stop("Reduced dimensions slot is NULL. Use runUmap().")
  }

  # custom hdbscan function to pass to fpc::clusterboot
  hdbscanCBI <-
    function(data, minPts, diss = inherits(data, "dist"), ...) {
      if (diss)
        c1 <- dbscan::hdbscan(data, minPts, method = "dist", ...)
      else
        c1 <- dbscan::hdbscan(data, minPts, ...)

      # check in case all of the cells are classified as outliers.
      if(length(unique(c1$cluster))==1){
        c1$cluster <- c1$cluster+1
      }

      tmp_data_df <- as.data.frame(data)
      tmp_data_df$hdb <- c1$cluster
      tmp_data_df$cell <- rownames(tmp_data_df)

      if (identical(rownames(tmp_data_df), rownames(data))) {
        c1$cluster <- tmp_data_df$hdb
        partition <- c1$cluster
      } else
        stop("order of dataframes is not identical.")

      #  plot(c1, data)

      cl <- list()
      nccl <- max(partition)
      partition[partition == 0] <- nccl + 1
      nc <- max(partition)
      #  print(nc)
      #  print(sc1)
      for (i in 1:nc)
        cl[[i]] <- partition == i
      out <-
        list(
          result = c1,
          nc = nc,
          nccl = nccl,
          clusterlist = cl,
          partition = partition,
          clustermethod = "hdbscan"
        )
      out
    }

  message(cat("Calculating jaccard similarity for k range:", k_range))

  jaccard <- BiocParallel::bplapply(k_range, function(i) {
    df_clusterboot <-
      .quiet(
        fpc::clusterboot(
          SingleCellExperiment::reducedDim(scCNA, "umap"),
          B = B,
          clustermethod = hdbscanCBI,
          seed = seed,
          minPts = i
        )
      )

    # -1 since the first partition is the outlier partition
    # which will be subclone 'c0'
    n_subclones <- sort(unique(df_clusterboot$partition)) - 1

    df_result <- data.frame(k = i,
                            subclones = paste0('c', n_subclones),
                            bootmean = df_clusterboot$bootmean)

  }, BPPARAM = BPPARAM)

  # Obtaining values from the list and binding to a dataframe
  names(jaccard) <- k_range
  bootmean_df <- do.call(rbind, jaccard)
  jc_df <- bootmean_df %>%
    dplyr::group_by(k) %>%
    dplyr::summarise(mean = mean(bootmean),
              median = median(bootmean))


  # Using selected summarising function
  if (fun == 'mean') {

    jc_df_opt <- jc_df %>%
      dplyr::filter(mean == max(jc_df$mean))

    selected_k_jac_value <- jc_df_opt$mean

  }

  if (fun == 'median') {

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
    "with", fun, "jaccard similarity of:",
    round(selected_k_jac_value, 3)
  ))

  # adding values to metadata
  S4Vectors::metadata(scCNA)$suggestedK_df <- bootmean_df
  S4Vectors::metadata(scCNA)$suggestedK <- jc_df_opt$k

  return(scCNA)

}
