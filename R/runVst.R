#' Variance Stabilizing Transformation
#'
#' Performs variance stabilization by using the Freeman-Tukey transformation
#'
#' @param scCNA The scCNA object
#' @param transformation Character. Transformation to be performed, available options are 'logratio' or 'freemantukey'
#'
#' @return A slot into the scCNA object containing the variance stabilized matrix. Can be accessed with \code{assay(scCNA, 'vst')}
#' @importFrom SummarizedExperiment assay
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

  if (transformation == 'logratio') {
    counts_df_ft <- purrr::map_dfc(varbin_counts_df, function(x) log2(x +1e-3))
  }

  counts_df_ft <- as.data.frame(counts_df_ft)

  assay(scCNA, transformation) <- counts_df_ft

  return(scCNA)

}
