#' Calculates the integer copy number profile for each single-cell
#'
#' @param scCNA The scCNA object.
#' @param assay String with the name of the assay to pull data from to calculate integers.
#' @param method Method used to scale the ratio values to integer.
#' @param ploidy_value If method of choice is 'fixed' a ploidy value should be provided.
#' @param name String specifying the name to be used to store the result in the
#' reducedDims of the output.
#'
#' @details
#' \itemize{
#' \item{fixed:} When method argument is set to 'fixed' copykit extracts the segment means
#' from the scCNA object and multiplies those means by the value provided in the
#' argument ploidy_value.
#' }
#'
#' @return The scCNA object with an assay slot named 'integer' that contains
#' a data frame with cells as columns and integerized segments as rows.
#' @export
#'
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay colData
#'
#' @examples
calcInteger <- function(scCNA,
                        assay = 'segment_ratios',
                        method = 'fixed',
                        ploidy_value = NULL,
                        name = 'integer') {

  seg_ratios_df <- SummarizedExperiment::assay(scCNA, assay)

  if (!is.null(ploidy_value)) {

    if (method == 'fixed') {

      if (is.null(ploidy_value) && !is.numeric(ploidy_value)) {
        stop("Method fixed requires a numeric value for ploidy_value argument.")
      }

      message(paste("Scaling ratio values by ploidy value",
                    ploidy_value))

      # ploidy values are added to colData information
      SummarizedExperiment::colData(scCNA)$ploidy <- ploidy_value

      # saving ploidy scaling method
      S4Vectors::metadata(scCNA)$ploidy_method <- 'fixed'

    }

  }

    # check to guarantee multiplication
    if (!identical(names(seg_ratios_df), colData(scCNA)$sample)) {
      stop("Order of cells in segment_ratios slot and colData(scCNA) is not identical.")
    }

    # obtain the matrix of integer values by multiplying the seg ratios
    # by the diagonal of the ploidy colData vector
    int_values <- round(as.matrix(seg_ratios_df) %*% diag(colData(scCNA)$ploidy)) %>%
      as.data.frame()

    # recovering names
    names(int_values) <- names(seg_ratios_df)

    SummarizedExperiment::assay(scCNA, assay) <- int_values

    return(scCNA)

}

