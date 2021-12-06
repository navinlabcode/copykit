#' calcInteger()
#'
#' Calculates the integer copy number profile for each single-cell
#'
#' @param scCNA The CopyKit object.
#' @param assay String with the name of the assay to pull data from to calculate
#' integers.
#' @param method Method used to scale the ratio values to integer.
#' @param ploidy_value If method of choice is 'fixed' a ploidy value should be
#' provided.
#' @param name String specifying the name to be used to store the result in the
#' reducedDims of the output.
#' @param penalty An integer passed on to scquantum::ploidy.inference()
#' penalty argument
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details
#' \itemize{
#' \item{fixed:} When method argument is set to 'fixed' copykit extracts the
#' segment means from the scCNA object and multiplies those means by the value
#' provided in the argument ploidy_value.
#' }
#'
#' @return The scCNA object with an assay slot named 'integer' that contains
#' a data frame with cells as columns and integerized segments as rows.
#' @export
#'
#' @importFrom S4Vectors metadata
#' @importFrom scquantum ploidy.inference
#' @importFrom SummarizedExperiment assay colData rowRanges
#'
#' @examples
#' copykit_obj <- mock_bincounts(ncells_diploid = 0)
#' copykit_obj <- calcInteger(copykit_obj, method = "scquantum")

calcInteger <- function(scCNA,
                        assay = "bincounts",
                        method = "fixed",
                        ploidy_value = NULL,
                        name = "integer",
                        penalty = 25,
                        BPPARAM = bpparam()) {

    df <- SummarizedExperiment::assay(scCNA, assay)

    if (!is.null(ploidy_value)) {
        if (method == "fixed") {
            if (is.null(ploidy_value) && !is.numeric(ploidy_value)) {
                stop("Method fixed requires a numeric value for ploidy_value.")
            }

            message(
                "Scaling ratio values by ploidy value ",
                ploidy_value
            )

            # ploidy values are added to colData information
            SummarizedExperiment::colData(scCNA)$ploidy <- ploidy_value

            # saving ploidy scaling method
            S4Vectors::metadata(scCNA)$ploidy_method <- "fixed"
        }
    }

    if (method == "scquantum") {
        rg <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))

        sc_quants <-
            BiocParallel::bplapply(
                assay(scCNA, assay),
                scquantum::ploidy.inference,
                chrom = rg$seqnames,
                start = rg$start,
                end = rg$end,
                penalty = penalty,
                BPPARAM = BPPARAM
            )

        sc_ploidies <- sapply(sc_quants, function(x) x$ploidy)
        sc_confidence <- sapply(sc_quants, function(x) x$confidence_ratio)

        SummarizedExperiment::colData(scCNA)$ploidy <- sc_ploidies
        SummarizedExperiment::colData(scCNA)$ploidy_confidence <- sc_confidence
    }

    # check to guarantee multiplication
    if (!identical(names(df), colData(scCNA)$sample)) {
        stop("Order of cells in segment_ratios and colData() is not identical.")
    }

    # obtain the matrix of integer values by multiplying the seg ratios
    # by the diagonal of the ploidy colData vector
    int_values <-
        round(as.matrix(df) %*% diag(colData(scCNA)$ploidy)) %>%
        as.data.frame()

    # recovering names
    names(int_values) <- names(df)

    SummarizedExperiment::assay(scCNA, name) <- int_values

    return(scCNA)
}
