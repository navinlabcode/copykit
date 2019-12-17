#' Filter noise cells
#'
#' filterCells uses a k-nearest-neighbor approach to remove cells with random CNA profiles, largely due to noise data. It calculates a correlation matrix and sets a resolution below which non neighbors will be classified as noise cells.
#'
#' @author Hua-Jun Wu
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param k K-nearest-neighbor, defaults to 5.
#' @param resolution Set's how strict the correlation cut off will be. Defaults to 0.8.
#' @param n_threads Number of parallel threads to calculate distances with \code{amap::Dist()}. Defaults to 1
#'
#' @return Adds a filtered cells label to the scCNA metadata. Cells that pass the filtering criteria receive the label "kept", whereas cells that do not pass the filtering criteria receive the label "removed".
#' @return Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @export
#'
#' @examples
#'
#'

filterCells <- function(scCNA,
                        k = 5,
                        resolution = 0.8,
                        n_threads = 1) {
  k = k
  seg <- copykit::segment_ratios(scCNA)

  message("Calculating correlation matrix.")
  dst = cor(seg)
  dst_knn_df = apply(as.matrix(dst), 1, function(x) {
    mean(sort(x, decreasing = T)[2:(k + 1)])
  }) %>%
    tibble::enframe(name = "sample",
                    value = "cor")

  dst_knn_df <- dst_knn_df %>%
    dplyr::mutate(filtered = dplyr::case_when(cor >= resolution ~ "kept",
                                              cor < resolution ~ "removed"))

  message(
    "Adding information to metadata. Access with SummarizedExperiment::colData(scCNA)."
  )
  if (identical(SummarizedExperiment::colData(scCNA)$sample,
                dst_knn_df$sample)) {
    SummarizedExperiment::colData(scCNA)$filtered <- dst_knn_df$filtered
  } else
    stop("Sample names do not match metadata sample info. Check colData(scCNA).")

  message("Plotting heatmap.")

  if (nrow(dst_knn_df) > 500) {
    message(
      paste(
        "Your dataset has:",
        nrow(dst_knn_df),
        "Cells. Plotting heatmap may take a long time with large number of cells. Set number of threads with n_threads for parallel processing if possible to speed up."
      )
    )
  }

  # filtered annotation
  filter_anno <-
    ComplexHeatmap::rowAnnotation(filtered = dplyr::pull(dst_knn_df, filtered),
                                  col = list(filtered = c(
                                    "kept" = "green2",
                                    "removed" = "firebrick3"
                                  )))

  # Heatmap
  ht <- ComplexHeatmap::Heatmap(
    t(seg),
    cluster_rows = function(x) {
      fastcluster::hclust(amap::Dist(x, method = "manhattan", nbproc = n_threads),
                          method = "ward.D2")
    },
    cluster_columns = FALSE,
    use_raster = TRUE,
    border = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_split = dplyr::pull(dst_knn_df, filtered),
    left_annotation = filter_anno
  )

  print(ht)

  message("Done.")
  return(scCNA)

}
