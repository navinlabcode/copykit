#' Calculates the ratios from a matrix of counts
#'
#' @param scCNA The scCNA object
#' @param assay String with the name of the assay to pull data from to calculate
#' the ratios.
#' @param fun A string indicating the summarizing function to be used.
#'
#' @details Calculates a sample-wise normalization of the selected assay by the
#' mean bin counts returns ratios where a value of 1 corresponds to the neutral
#'copy number state of the sample
#'
#' @return A ratio matrix within the slot assay(scCNA, 'ratios')
#' can be accessed with \code{ratios}.
#' @export
#'
#' @importFrom SummarizedExperiment assay
#'
#' @examples
calcRatios <- function(scCNA,
                       assay = c("ft", "bin_counts"),
                       fun = c("mean", "median")) {

  if (assay %!in% c("ft", "bin_counts")) {
    stop("Assay must be either 'ft' or 'bin_counts'")
  }

  assay <- match.arg(assay)
  fun <- match.arg(fun)

  counts <- SummarizedExperiment::assay(scCNA, assay)

  ratios_df <- sweep(counts, 2, apply(counts, 2, fun), '/')

  SummarizedExperiment::assay(scCNA, "ratios") <- ratios_df

  return(scCNA)

}
