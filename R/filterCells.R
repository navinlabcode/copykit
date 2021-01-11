#' Filter noise cells
#'
#' filterCells uses a k-nearest-neighbor approach to remove cells
#' with random CNA profiles, largely due to noise data.
#' It calculates a correlation matrix and sets a resolution below which non neighbors
#'  will be classified as noise cells.
#'
#' @author Hua-Jun Wu
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param k K-nearest-neighbor, defaults to 5.
#' @param resolution Set's how strict the correlation cut off will be. Defaults to 0.8.
#' @param n_threads Number of parallel threads to calculate distances with \code{amap::Dist()}. Defaults to 1/4 of the cores available in your system.
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
                        n_threads = parallel::detectCores() / 4) {
  # checks
  if (n_threads < 1) {
    n_threads <- 1
  }

  if (!is.numeric(resolution)) {
    stop("Resolution needs to be a number between 0 and 1")
  }

  if (resolution < 0 || resolution > 1) {
    stop("Resolution needs to be a number between 0 and 1")
  }

  seg <- copykit::segment_ratios(scCNA)

  message("Calculating correlation matrix.")
  dst <- cor(seg)
  dst_knn_df <- apply(as.matrix(dst), 1, function(x) {
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
    SummarizedExperiment::colData(scCNA)$filter_corr_value <- dst_knn_df$cor
    SummarizedExperiment::colData(scCNA)$filtered <- dst_knn_df$filtered

  } else
    stop("Sample names do not match metadata sample info. Check colData(scCNA).")

  message("Done.")
  return(scCNA)

}
