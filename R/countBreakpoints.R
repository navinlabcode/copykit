#' countBreakpoints
#'
#' Considers changes in the segment ratios as breakpoints.
#' Counts the breakpoints for each chromosome arm separately.
#'
#' @param scCNA
#'
#' @return The scCNA object with a column of breakpoint counts added to the colData.
#' @export
#'
#' @keywords internal
#'
#' @examples
.countBreakpoints <- function(scCNA) {

  rg_chr <- SummarizedExperiment::rowRanges(scCNA) %>%
    as.data.frame() %>%
    dplyr::mutate(chrarm = paste0(seqnames, arm)) %>%
    dplyr::select(chrarm)

  dat_seg_cp <- segment_ratios(scCNA)

  # split by chrom
  message("Counting breakpoints.")
  dat_seg_split <- split(dat_seg_cp, pull(rg_chr, chrarm))

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

  return(scCNA)
}
