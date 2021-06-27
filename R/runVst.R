#' Variance Stabilizing Transformation
#'
#' Performs variance stabilization by using the Freeman-Tukey transformation
#'
#' @param scCNA The scCNA object
#' @param transformation Character. Transformation to be performed, available options are 'log' or 'ft'
#'
#' @details \code{runVst} performs variance stabilization to reduce the overdispersion
#' from the negative binomial distribution nature of the bin counts and reduce
#' technical bias. The argument \code{vst} controls the choice of the transformation
#' allowing either the Freeman-Tukey transformation by using the option 'ft' (recommended)
#' or a logarithmic transformation with the option 'log'. Using a 'log' transformation
#' may result in long segmentation times for a few cells with large breakpoint counts.
#'
#' @return A slot into the scCNA object containing the variance stabilized matrix.
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @importFrom purrr map_dfc
#' @export
#'
#' @examples
runVst <- function(scCNA,
                   transformation = 'ft') {

  varbin_counts_df <- copykit::bin_counts(scCNA)

  if (transformation == 'ft') {
    counts_df_ft <- purrr::map_dfc(varbin_counts_df, function(x) sqrt(x) + sqrt(x+1))
  }

  if (transformation == 'log') {
    counts_df_ft <- purrr::map_dfc(varbin_counts_df, function(x) log2(x+1e-3))
  }

  counts_df_ft <- as.data.frame(counts_df_ft)

  S4Vectors::metadata(scCNA)$vst <- transformation
  SummarizedExperiment::assay(scCNA, transformation) <- counts_df_ft

  return(scCNA)

}
