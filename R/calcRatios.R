#' Calculates the ratios from a matrix of counts
#'
#' @param scCNA The scCNA object
#' @param assay The assay that will be used
#' @param fun Character. Function used to calculate the ratios, defaults to "median"
#'
#' @return A ratio matrix within the slot \code{assay(scCNA, 'ratios')} can be accessed with \code{copykit::ratios(scCNA)}.
#' @export
#'
#' @importFrom SummarizedExperiment assay
#'
#' @examples
calcRatios <- function(scCNA,
                       assay = "ft",
                       fun = "median") {

  if (assay %!in% c("ft", "bincounts")) {
    stop("Assay must be either 'ft' or 'bincounts'" )
  }

  counts <- assay(scCNA, assay)

  ratios_df <- sweep(counts, 2, apply(counts, 2, fun), '/')

  assay(scCNA, "ratios") <- ratios_df

  return(scCNA)

}
