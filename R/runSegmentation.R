#' Run Segmentation
#'
#' Runs a segmentation algorithm using the ratio data.
#'
#' @param scCNA The scCNA object
#' @param method Character. Segmentation method of choice.
#'
#' @return The segment profile for all cells inside the scCNA object. Can be retrieved with \link{\code{copykit::segment_ratios()}}
#' @export
#'
#' @examples
runSegmentation <- function(scCNA,
                            method = "CBS") {
  ratios_df <- copykit::ratios(scCNA)

  if (method == "CBS") {
    CBS_seg <- parallel::mclapply(ratios_df, function(x) {
      CNA_object <-
        DNAcopy::CNA(
          log(x + 1e-3, base = 2),
          as.numeric(stringr::str_remove(hg38_rg$chr, "chr")),
          hg38_rg$start,
          data.type = "logratio",
          sampleid = names(x)
        )
      smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
      segment_smoothed_CNA_object <-
        DNAcopy::segment(
          smoothed_CNA_object,
          alpha = 0.01,
          min.width = 5,
          undo.splits = "prune",
          undo.prune = 0.05
        )
      short_cbs <- segment_smoothed_CNA_object[[2]]
      log_seg_mean_LOWESS <-
        rep(short_cbs$seg.mean, short_cbs$num.mark)
      merge_obj <-
        MergeLevels(smoothed_CNA_object[, 3], log_seg_mean_LOWESS)$vecMerged
      merge_ratio <- 2 ^ merge_obj

    }, mc.cores = n_threads)

    cbs_seg_df <- bind_cols(CBS_seg) %>%
      as.data.frame()

    copykit::segment_ratios(scCNA) <- cbs_seg_df

  }


}

