#' Calculates the integer copy number profile for each single-cell
#'
#' @param scCNA The scCNA object.
#' @param method Method used to calculate integer profiles. Options: 'fixed'.
#' @param ploidy_value If method of choice is 'fixed' a ploidy value should be provided
#'
#' @details
#' \item{fixed:} When method argument is set to 'fixed' copykit extracts the segment means
#' from the scCNA object and multiplies those means by the value provided in the
#' argument ploidy_value.
#'
#' @return The scCNA object with an assay slot named 'integer' that contains
#' a data frame with cells as columns and integerized segments as rows.
#' @export
#'
#' @examples
calcInteger <- function(scCNA,
                        method,
                        ploidy_value = NULL) {
  seg_ratios_df <- copykit::segment_ratios(scCNA)

  if (method == "fixed") {
    if (is.null(ploidy_value) && !is.numeric(ploidy_value)) {
      stop("Method fixed requires a numeric value for ploidy_value argument.")
    }

    message(paste("Scaling ratio values by ploidy value",
                  ploidy_value))

    int_values <- round(seg_ratios_df * ploidy_value)

    assay(scCNA, 'integer') <- int_values

    return(scCNA)

  }

}

