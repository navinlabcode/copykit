#' Variance Stabilizing Transformation
#'
#' Performs variance stabilization transformation of the bin counts
#'
#' @param scCNA The scCNA object
#' @param transformation A character indicating the variance stabilization
#' transformation to be performed. See \link{runVst} details.
#'
#' @details \code{runVst} performs variance stabilization to reduce the overdispersion
#' from the negative binomial distribution nature of the bin counts and reduce
#' technical bias. The argument \code{vst} controls the choice of the transformation
#' allowing either the Freeman-Tukey transformation by using the option 'ft' (recommended)
#' or a logarithmic transformation with the option 'log'. Using a 'log' transformation
#' may result in long segmentation times for a few cells with large breakpoint counts.
#'
#' @references
#' Freeman, M. F.; Tukey, J. W. (1950), "Transformations related to the angular
#' and the square root", The Annals of Mathematical Statistics,
#' 21 (4), pp. 607–611, doi:10.1214/aoms/1177729756, JSTOR 2236611
#'
#' @return A slot into the scCNA object containing the variance stabilized matrix.
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors metadata
#' @importFrom purrr map_dfc
#' @export
#'
#' @examples
runVst <- function(scCNA,
                   transformation = c('ft','log')) {

  transformation <- match.arg(transformation)

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
