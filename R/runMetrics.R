#' Run metrics
#'
#' Calculates the RMSE and the breakpoint counts for each cell.
#' Stores it in the metadata that can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' Results can be visualized with \code{plotMetrics()}
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @return Adds a the metrics to the scCNA metadata. Those metrics can be used for subsetting the data if desired
#' @return Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @export
#' @import ggplot2
#'
#' @examples
#'
#'
runMetrics <- function(scCNA,
                       BPPARAM = bpparam()) {

  ###################
  # Retrieving data

  dat_seg <- copykit::segment_ratios(scCNA)
  dat_rat <- copykit::ratios(scCNA)
  rg <- SummarizedExperiment::rowRanges(scCNA)

  ####################
  # RMSE
  message("Calculating RMSE")

  rmses <- BiocParallel::bplapply(1:ncol(dat_seg), function(i) {
    actual <- dat_seg[, i]
    predicted <- dat_rat[, i]
    Metrics::rmse(actual, predicted)
  },
  BPPARAM = BPPARAM)

  rmses <- unlist(rmses)

  SummarizedExperiment::colData(scCNA)$rmse <- rmses

  ####################
  # Number of breakpoints
  # Performed for every chromosome

  rg_chr <- rg %>%
    as.data.frame() %>%
    dplyr::select(seqnames) %>%
    dplyr::mutate(seqnames = as.character(seqnames))

  dat_seg_cp <- dat_seg

  # split by chrom
  message("Counting breakpoints.")
  dat_seg_split <- split(dat_seg_cp, pull(rg_chr, seqnames))

  brkpt_by_chrom <-
    lapply(dat_seg_split, function(x) {
      purrr::map_dfc(x, function(i) {
        length(rle(i)$values) - 1
      }) %>%
        unlist()
    })

  brkpt_by_chrom_df <- dplyr::bind_rows(brkpt_by_chrom) %>%
    t() %>%
    as.data.frame()

  brkpt_by_chrom_l <- brkpt_by_chrom_df %>%
    tibble::rownames_to_column(var = "sample") %>%
    tidyr::gather(key = "chr",
                  value = "brkpts",-sample)

  brkpt_by_sample_cnt <- brkpt_by_chrom_l %>%
    dplyr::group_by(sample) %>%
    dplyr::tally(brkpts)

  # making sure order is correct
  brkpt_by_sample_cnt <-
    brkpt_by_sample_cnt[match(brkpt_by_sample_cnt$sample,
                              SummarizedExperiment::colData(scCNA)$sample), ]

  SummarizedExperiment::colData(scCNA)$breakpoint_count <-
    brkpt_by_sample_cnt$n

  message("Done.")

  return(scCNA)

}
