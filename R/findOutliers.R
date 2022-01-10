#' findOutliers()
#'
#' Uses a nearest neighbor approach to find noise copy number profiles within
#' the segment means.
#'
#' @author Hua-Jun Wu
#' @author Darlan Conterno Minussi
#' @author Junke Wang
#'
#' @param scCNA CopyKit object.
#' @param assay String with the name of the assay to pull data.
#' @param k A numeric scalar with the number k-nearest-neighbor cells to
#' calculate the mean correlation
#' @param resolution A numeric scalar that set's how strict the
#' correlation cut off will be.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details \code{findOutliers} To detect low-quality cells, CopyKit calculates
#'  the Pearson correlation matrix of all samples from the segment ratio means.
#'  Next, we calculate a sample-wise mean of the correlation between a sample
#'  and its k-nearest-neighbors. Samples in which the correlation value is lower
#'  than the defined threshold are classified as low-quality cells.
#'
#' @return Adds a column 'outlier' to
#' \code{\link[SummarizedExperiment]{colData}}. Cells that pass the filtering
#' criteria receive the label "kept", whereas cells that do not pass the
#' filtering criteria receive the label "removed".
#'
#' @importFrom stats cor sd
#'
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example()
#' copykit_obj <- findNormalCells(copykit_obj)
#' copykit_obj <- copykit_obj[, colData(copykit_obj)$is_normal == "FALSE"]
#' copykit_obj <- findOutliers(copykit_obj)
findOutliers <- function(scCNA,
    assay = "segment_ratios",
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

    dst <- parCor(seg, BPPARAM = BPPARAM)

    dst_knn_df <- apply(as.matrix(dst), 1, function(x) {
        mean(sort(x, decreasing = TRUE)[2:(k + 1)])
    }) %>%
        tibble::enframe(
            name = "sample",
            value = "cor"
        )

    dst_knn_df <- dst_knn_df %>%
        dplyr::mutate(outlier = dplyr::case_when(
            cor >= resolution ~ "FALSE",
            cor < resolution ~ "TRUE"
        ))

    n_filtered <- table(dst_knn_df$outlier)["TRUE"]
    message("Marked ", n_filtered, " cells as outliers.")

    message(
        "Adding information to metadata. Access with colData(scCNA)."
    )
    if (identical(
        SummarizedExperiment::colData(scCNA)$sample,
        dst_knn_df$sample
    )) {
        SummarizedExperiment::colData(scCNA)$cell_corr_value <-
            round(dst_knn_df$cor, 3)
        SummarizedExperiment::colData(scCNA)$outlier <-
            dst_knn_df$outlier
    } else {
        stop("Sample names do not match colData info. Check colData(scCNA).")
    }

    message("Done.")
    return(scCNA)
}
