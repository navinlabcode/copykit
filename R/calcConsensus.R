#' Calculate a consensus matrix of segment means based on \code{colData}
#'
#' @param scCNA The CopyKit object.
#' @param assay String with the name of the assay to pull data from to calculate
#' the consensus matrix.
#' @param consensus_by A string with the column from colData that will be used
#'  to isolate the cells by factor and calculate the consensus.
#' @param fun A string indicating the summarizing function to be used.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details Consensus profiles are calculated by averaging or taking the median
#'  of the ith segment mean of all single cells assigned to the same element of
#'  \link{colData},
#'
#' @return A consensus matrix stored in the consensus slot of the CopyKit object
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- findClusters(copykit_obj)
#' copykit_obj <- calcConsensus(copykit_obj)
calcConsensus <- function(scCNA,
                          assay = "segment_ratios",
                          consensus_by = "subclones",
                          fun = c("median", "mean"),
                          BPPARAM = bpparam()) {
    fun <- match.arg(fun)

    if (consensus_by == "subclones" &
        is.null(SummarizedExperiment::colData(scCNA)$subclones)) {
        stop("Calculating consensus requires clusters. use findClusters(scCNA)")
    }

    if (consensus_by %!in% names(SummarizedExperiment::colData(scCNA))) {
        stop("consensus_by must be an element of colData(scCNA)")
    }

    if (length(consensus_by) != 1) {
        stop("consensus_by argument must have length == 1")
    }

    if (is.null(consensus_by)) {
        stop("Please provide information to consensus_by argument.")
    }

    consensus_info <-
        as.data.frame(SummarizedExperiment::colData(scCNA)) %>%
        dplyr::select(!!consensus_by) %>%
        droplevels()

    seg_data <- as.data.frame(t(SummarizedExperiment::assay(scCNA, assay)))

    # sanity check
    if (!identical(rownames(consensus_info), rownames(seg_data))) {
        stop("Order of elements in colData and segment_ratios must be identical.")
    }

    ## reading list with clusters
    long_list <- split(seg_data, consensus_info)

    consensus_list <-
        BiocParallel::bplapply(long_list, function(x) {
            apply(x, 2, fun)
        }, BPPARAM = BPPARAM)

    cs_df <- as.data.frame(t(do.call(rbind, consensus_list)))

    if (assay == "integer") {
        cs_df <- round(cs_df)
    }

    names(cs_df) <- names(consensus_list)

    # This hidden attribute will allow plotHeatmap to figure it out which
    # argument was used in 'consensus_by'
    attr(cs_df, "consensus_by") <- consensus_by

    # This hidden attribute will allow plotHeatmap to figure it out which
    # argument was used in 'assay'
    attr(cs_df, "consensus_assay") <- assay

    consensus(scCNA) <- cs_df

    return(scCNA)
}
