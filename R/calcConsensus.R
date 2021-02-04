#' Calculate a consensus matrix of segment ratios based on metadata
#'
#' @param scCNA The scCNA object.
#' @param consensus_by The column from metadata that will be used to isolate the cells by factor and calculate the consensus.
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `parallel::mclapply`. As default it uses 1/4 of the detected cores available.
#'
#' @return
#' @export
#'
#' @examples
calcConsensus <- function(scCNA,
                          consensus_by = "subclones",
                          n_threads = parallel::detectCores() / 4) {

  if (consensus_by == 'subclones' & is.null(SummarizedExperiment::colData(scCNA)$subclones)) {
    stop("Calculating consensus requires cluster information. use findClusters(scCNA)")
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

  consensus_info <- as.data.frame(SummarizedExperiment::colData(scCNA)) %>%
    dplyr::select(!!consensus_by)

  seg_data <- as.data.frame(t(segment_ratios(scCNA)))

  #sanity check
  if (!identical(rownames(consensus_info), rownames(seg_data))) {
    stop("Order of elements in metadata and in segment_ratios information must be identical.")
  }

  ## reading list with clusters
  long_list <- split(seg_data, consensus_info)

  consensus_list <-
    parallel::mclapply(long_list, function(x) {
      apply(x, 2, median)
    }, mc.cores = n_threads)

  cs_df <- as.data.frame(t(do.call(rbind, consensus_list)))

  names(cs_df) <- names(consensus_list)

  # This hidden attribute will allow plotHeatmap to figure it out which
  # argument was used in 'consensus_by'
  attr(cs_df, "consensus_by") <- consensus_by

  consensus(scCNA) <- cs_df

  return(scCNA)

}
