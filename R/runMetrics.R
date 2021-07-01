#' Run metrics
#'
#' Calculates the overdispersion and the breakpoint counts for each cell.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details Adds the metrics to the scCNA \code{\link[SummarizedExperiment]{colData}}.
#'  Those metrics can be used for subsetting the data if desired.
#'  results can be visualized with \code{\link{plotMetrics}}.
#'
#' @return Adds columns 'overdispersion' and 'breakpoint_count' to
#'  \code{\link[SummarizedExperiment]{colData}}.
#'
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

  dat_seg <- segment_ratios(scCNA)
  dat_rat <- ratios(scCNA)
  dat_bin <- bin_counts(scCNA)
  rg <- SummarizedExperiment::rowRanges(scCNA)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:01 2021
  # overdispersion ----
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:10 2021

  message("Calculating overdispersion.")

  overdisp <- BiocParallel::bplapply(dat_bin,
                                  overdispersion,
                                  BPPARAM = BPPARAM)

  overdisp <- unlist(overdisp)

  SummarizedExperiment::colData(scCNA)$overdispersion <- overdisp

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:25 2021
  # Breakpoint count
  # Performed for every chromosome
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:22:34 2021

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
