#' Run Segmentation
#'
#' Runs a segmentation algorithm using the ratio data.
#'
#' @param scCNA The scCNA object
#' @param method Character. Segmentation method of choice.
#' @param genome Character. Genome assembly to be used, current accepted "hg19" or "hg38".
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `parallel::mclapply`. As default it uses 1/4 of the detected cores available.
#'
#' @return The segment profile for all cells inside the scCNA object. Can be retrieved with \code{copykit::segment_ratios()}
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @importMethodsFrom SummarizedExperiment
#' @export
#'
#' @examples
runSegmentation <- function(scCNA,
                            method = "CBS",
                            genome = "hg38",
                            n_threads = parallel::detectCores() / 4) {

  # checks

  if (!is.numeric(n_threads)) {
    stop("n_threads argument must be numeric.")
  }

  if (n_threads > parallel::detectCores()) {
    stop(paste(
      "n_threads argument must be smaller than",
      parallel::detectCores()
    ))
  }

  if (n_threads < 1) {
    n_threads <- 1
  }

  if (genome %!in% c("hg19", "hg38")) {
    stop("Genome assembly must be 'hg19' or 'hg38'")
  }

  # genome assembly
  # Reading hg38 VarBin ranges
  if (genome == "hg38") {

    #match for chrY presence
    chr_sccna <- as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
    hg38_rg <- hg38_rg[which(hg38_rg$chr %in% chr_sccna),]

    hg38_rg <- hg38_rg %>%
      mutate(chr = str_replace(chr, "X", "23"),
             chr = str_replace(chr, "Y", "24"))

    chr_info <-  as.numeric(stringr::str_remove(hg38_rg$chr, "chr"))

    ref <- hg38_rg

  }

  # reading hg19 varbin ranges
  if (genome == "hg19") {

    #match for chrY presence
    chr_sccna <- as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
    hg19_rg <- hg19_rg[which(hg19_rg$chr %in% chr_sccna),]

    hg19_rg <- hg19_rg %>%
      mutate(chr = str_replace(chr, "X", "23"),
             chr = str_replace(chr, "Y", "24"))

    chr_info <-  as.numeric(stringr::str_remove(hg19_rg$chr, "chr"))

    ref <- hg19_rg

  }

  ratios_df <- copykit::ratios(scCNA)

  if (method == "CBS") {
    CBS_seg <- parallel::mclapply(ratios_df, function(x) {
      CNA_object <-
        DNAcopy::CNA(
          log2(x+1e-3),
          chr_info,
          ref$start,
          data.type = "logratio",
          sampleid = names(x)
        )
      smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
      segment_smoothed_CNA_object <-
        DNAcopy::segment(
          smoothed_CNA_object,
          alpha = 0.01,
          min.width = 5,
          undo.splits = "sdundo",
          undo.SD = 1,
          # undo.prune = 0.05
        )
      short_cbs <- segment_smoothed_CNA_object[[2]]
      log_seg_mean_LOWESS <-
        rep(short_cbs$seg.mean, short_cbs$num.mark)
      merge_obj <-
        .MergeLevels(smoothed_CNA_object[, 3], log_seg_mean_LOWESS)$vecMerged
      merge_ratio <- 2^merge_obj

    }, mc.cores = n_threads)

    cbs_seg_df <- bind_cols(CBS_seg) %>%
      as.data.frame()

    SummarizedExperiment::assay(scCNA, 'segment_ratios') <- cbs_seg_df

    return(scCNA)

  }


}

