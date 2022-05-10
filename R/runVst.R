#' Variance Stabilizing Transformation
#'
#' Performs variance stabilization transformation of the bin counts
#'
#' @param scCNA The scCNA object
#' @param transformation A character indicating the variance stabilization
#' transformation to be performed. See \link{runVst} details.
#' @param assay A character indicating the assay slot to extract the bincounts
#' for variance stabilization
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
#' @export
#'
#' @examples
#' copykit_obj <- mock_bincounts(ncells = 10)
#' copykit_obj <- runVst(copykit_obj)
runVst <- function(scCNA,
                   transformation = c("ft", "log"),
                   assay = 'bincounts') {
    transformation <- match.arg(transformation)

    message(paste("Running variance stabilization transformation:",
                  transformation))

    # recovering assay
    varbin_counts_df <- assay(scCNA, assay)

    if (transformation == "ft") {
        counts_df_ft <- as.data.frame(apply(varbin_counts_df,
                                            2,
                                            function(x) sqrt(x) + sqrt(x + 1)))
    }

    if (transformation == "log") {
        varbin_counts_df[varbin_counts_df == 0] <- 1e-4
        counts_df_ft <- as.data.frame(apply(varbin_counts_df,
                                            2,
                                            function(x) log(x)))
    }

    counts_df_ft <- as.data.frame(counts_df_ft)

    S4Vectors::metadata(scCNA)$vst <- transformation
    SummarizedExperiment::assay(scCNA, transformation) <- counts_df_ft

    return(scCNA)
}
