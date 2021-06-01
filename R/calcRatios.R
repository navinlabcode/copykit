#' Calculates the ratios from a matrix of counts
#'
#' @param scCNA The scCNA object
#' @param assay The assay that will be used
#' @param fun Character. Function used to calculate the ratios.
#' Defaults to "median"
#'
#' @return A ratio matrix within the slot \code{assay(scCNA, 'ratios')}
#' can be accessed with \code{copykit::ratios(scCNA)}.
#' @export
#'
#' @importFrom SummarizedExperiment assay
#'
#' @examples
calcRatios <- function(scCNA,
                       assay = "ft",
                       fun = "mean") {
  if (assay %!in% c("ft", "bin_counts")) {
    stop("Assay must be either 'ft' or 'bin_counts'")
  }

  counts <- SummarizedExperiment::assay(scCNA, assay)

  ratios_df <- sweep(counts, 2, apply(counts, 2, fun), '/')

  SummarizedExperiment::assay(scCNA, "ratios") <- ratios_df

  return(scCNA)

}
