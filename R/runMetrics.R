#' Run metrics
#'
#' Calculates the RMSE and the breakpoint counts for each cell.
#' Stores it in the metadata that can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' Results can be visualized with \code{plotMetrics()}
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `parallel::mclapply`. As default it uses 1/4 of the detected cores available.
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
                       n_threads = parallel::detectCores() / 4) {
  # checks

  if (!is.numeric(n_threads)) {
    stop("n_threads argument must be numeric.")
  }

  # checks
  if (n_threads > parallel::detectCores()) {
    stop(paste(
      "n_threads argument must be smaller than ",
      parallel::detectCores()
    ))
  }

  if (n_threads < 1) {
    n_threads <- 1
  }

  ###################
  # Retrieving data

  dat_seg <- copykit::segment_ratios(scCNA)
  dat_rat <- copykit::ratios(scCNA)
  rg <- SummarizedExperiment::rowRanges(scCNA)

  ####################
  # RMSE
  message("Calculating RMSE")
  message(paste("Using", n_threads, "cores."))

  rmses <- parallel::mclapply(1:ncol(dat_seg), function(i) {
    actual <- dat_seg[, i]
    predicted <- dat_rat[, i]
    Metrics::rmse(actual, predicted)
  },
  mc.cores = n_threads)

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
    as.data.frame()
  rownames(brkpt_by_chrom_df) <- colnames(dat_seg_cp)

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
