#' Finds the Optimal K value to be used for subclone clustering
#'
#' @param scCNA  scCNA object.
#' @param k_range Range of values to be tested. Defaults to 7 to the sqrt of the number of cells
#' @param method Method which where the values will be tested. Only "hdbscan" available.
#' @param seed Seed (Defaults to 17).
#' @param B Number of bootstrapping. Defaults to 100. Higher values yield better results at a cost of performance
#' @param n_threads  Number of threads. Passed to `parallel::mclapply`. As default it uses 1/4 of the detected cores available.
#'
#' @return
#' @export
#'
#' @importFrom fpc clusterboot
#' @importFrom tibble enframe
#' @importFrom dbscan hdbscan
#' @importFrom SingleCellExperiment reducedDim
#'
#' @examples
findOptimalK <- function(scCNA,
                         k_range = 7:sqrt(ncol(segment_ratios(scCNA))),
                         method = "hdbscan",
                         seed = 17,
                         B = 100,
                         n_threads = parallel::detectCores() / 4) {

  # obtaining data from reducedDim slot
  if (is.null(SingleCellExperiment::reducedDim(scCNA))) {
    stop("Reduced dimensions slot is NULL. Use runUmap().")
  }


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

      # for hdb
      # adding the ones classified as outliers to the closest cluster possible according to euclidean distance
      dist_umap <-
        dist(tmp_data_df[, c(1:2)]) %>% as.matrix() %>% as.data.frame() %>%
        rownames_to_column("cell2") %>%
        tidyr::gather(key = "cell1",
               value = "dist",-cell2) %>%
        dplyr::filter(cell1 != cell2)

      dist_min <- dist_umap %>%
        right_join(tmp_data_df %>% dplyr::select(cell, hdb),
                   by = c("cell2" = "cell")) %>%
        dplyr::filter(hdb != 0) %>%
        dplyr::group_by(cell1) %>%
        dplyr::slice_min(dist) %>%
        ungroup()


      for (i in 1:nrow(tmp_data_df)) {
        if (tmp_data_df$hdb[i] == 0) {
          cellname <- rownames(tmp_data_df)[i]
          closest_cell <-
            dplyr::filter(dist_min, cell1 == rownames(tmp_data_df)[i])$cell2
          closest_cell_cluster <-
            dplyr::filter(tmp_data_df, cell == closest_cell)$hdb
          tmp_data_df$hdb[i] <- closest_cell_cluster

        }

      }

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

  mean_jaccard <- parallel::mclapply(k_range, function(i) {
    df_clusterboot <-
      fpc::clusterboot(
        SingleCellExperiment::reducedDim(scCNA, "umap"),
        B = B,
        clustermethod = hdbscanCBI,
        seed = seed,
        minPts = i
      )

    mean_jc_cl <- mean(df_clusterboot$bootmean)

  }, mc.cores = n_threads)

  names(mean_jaccard) <- k_range

  mean_jc_df <- tibble::enframe(unlist(mean_jaccard),
                                name = "k",
                                value = "mean_jaccard") %>%
    mutate(k = as.numeric(k))

  mean_jc_df <- mean_jc_df %>%
    filter(mean_jaccard == max(mean_jc_df$mean_jaccard))

  if (nrow(mean_jc_df) > 1) {
    mean_jc_df <- mean_jc_df %>%
      filter(k == max(k))
  }

  message(paste(
    "Optimal k =",
    mean_jc_df$k,
    "with mean jaccard similarity of:",
    round(mean_jc_df$mean_jaccard, 3)
  ))

  return(mean_jc_df)

}
