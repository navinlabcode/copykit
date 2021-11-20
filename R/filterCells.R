#' Filter noise cells
#'
#' Uses a nearest neighbor approach to find noise copy number profiles within the
#' segment means.
#'
#' @author Hua-Jun Wu
#' @author Darlan Conterno Minussi
#' @author Junke Wang
#'
#' @detail \code{filterCells} Calculates a correlation matrix across the segment
#' means among all cells and takes the mean of its k-nearest neighbors correlation.
#' A threshold (argument resolution) is used for the minimum acceptable mean
#' correlation among the cell and its neighbors. Values below the set resolution
#' will be classified as noise cells.
#'
#' @param scCNA scCNA object.
#' @param assay String with the name of the assay to pull data.
#' @param k A numeric scalar with the number k-nearest-neighbor cells to
#' calculate the mean correlation
#' @param resolution A numeric scalar that set's how strict the
#' correlation cut off will be.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#'should be parallelized.
#'
#' @return Adds a column 'filtered' to \code{\link[SummarizedExperiment]{colData}}
#' Cells that pass the filtering criteria receive the label "kept",
#' whereas cells that do not pass the filtering criteria
#' receive the label "removed".
#'
#' @importFrom stats cor
#'
#' @export
#'
#' @examples
#'
#'

filterCells <- function(scCNA,
                        assay = 'segment_ratios',
                        k = 5,
                        resolution = 0.9,
                        BPPARAM = BiocParallel::bpparam()) {

  if (!is.numeric(resolution)) {
    stop("Resolution needs to be a number between 0 and 1")
  }

  if (resolution < 0 || resolution > 1) {
    stop("Resolution needs to be a number between 0 and 1")
  }

  seg <- SummarizedExperiment::assay(scCNA, assay)

  message("Calculating correlation matrix.")

  # correction to avoid correlations calculations with standard deviation zero
  zero_sd_idx <- which(apply(seg, 2, sd) == 0)

  if (length(zero_sd_idx) >= 1) {
    seg[1, zero_sd_idx] <- seg[1, zero_sd_idx] + 1e-3
  }

  # calculating correlations

  dst <- parCor(seg, BPPARAM=BPPARAM)

  dst_knn_df <- apply(as.matrix(dst), 1, function(x) {
    mean(sort(x, decreasing = T)[2:(k + 1)])
  }) %>%
    tibble::enframe(name = "sample",
                    value = "cor")

  dst_knn_df <- dst_knn_df %>%
    dplyr::mutate(filtered = dplyr::case_when(cor >= resolution ~ "kept",
                                              cor < resolution ~ "removed"))

  n_filtered <- table(dst_knn_df$filtered)['removed']
  message(paste("Marked", n_filtered, "cells to be removed."))

  message(
    "Adding information to metadata. Access with colData(scCNA)."
  )
  if (identical(SummarizedExperiment::colData(scCNA)$sample,
                dst_knn_df$sample)) {
    SummarizedExperiment::colData(scCNA)$filter_corr_value <-
      round(dst_knn_df$cor, 3)
    SummarizedExperiment::colData(scCNA)$filtered <-
      dst_knn_df$filtered

  } else
    stop("Sample names do not match metadata sample info. Check colData(scCNA).")

  message("Done.")
  return(scCNA)

}
