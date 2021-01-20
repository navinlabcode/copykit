#' Calculates the ratios from a matrix of counts
#'
#' @param scCNA The scCNA object
#' @param assay The assay that will be used
#' @param fun Character. Function used to calculate the ratios, defaults to "median"
#'
#' @return A ratio matrix within the slot \link{\code{assay(scCNA, 'ratios')}} can be accessed with \link{\code{copykit::ratios(scCNA)}}.
#' @export
#'
#' @examples
calcRatios <- function(scCNA,
                       assay = "vst",
                       fun = "median") {

  counts <- assay(scCNA, assay)

  ratios_df <- sweep(counts, 2, apply(counts, 2, fun), '/')

  assay(scCNA, "ratios") <- ratios_df

  return(scCNA)

}
