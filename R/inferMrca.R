#' inferMrca
#'
#' From the consensus matrix it infers a Most Recent Common Ancestral (MRCA)
#' across all groups.
#'
#' @param scCNA  the scCNA object
#' @param value  A numeric value used to compare the profiles to infer the mrca,
#' usually equal to the mean segment ratio of cells (value = 1) or the average
#' copy number of the cells
#'
#' @details Calculates the MRCA by inferring, for every bin, the value across
#' all groups that is closest to the number supplied in the argument value.
#'
#' @return Returns a numeric vector added to the \code{\link[S4Vectors]{metadata}}
#' of the scCNA object named `inferred_mrca`
#' @export
#'
#' @importFrom SummarizedExperiment seqnames
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()[,1:300]
#' copykit_obj <- findClusters(copykit_obj)
#' copykit_obj <- calcConsensus(copykit_obj)
#' copykit_obj <- inferMrca(copykit_obj)
inferMrca <- function(scCNA,
                      value = 1) {
    if (nrow(consensus(scCNA)) == 0) {
        stop("Consensus slot is empty. run calcConsensus().")
    }

    consensus_df <- as.data.frame(t(consensus(scCNA)))

    anc_profile <- apply(
        consensus_df,
        2,
        function(x) x[which.min(abs(x - value))]
    )

    metadata(scCNA)$inferred_mrca <- anc_profile

    return(scCNA)
}
