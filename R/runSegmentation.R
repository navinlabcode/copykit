#' Run Segmentation
#'
#' Runs a segmentation algorithm using the ratio data.
#'
#' @param scCNA The scCNA object
#' @param method Character. Segmentation method of choice.
#' @param genome Character. Genome assembly to be used, current accepted "hg19" or "hg38".
#' @param seed Numeric. Set seed for CBS segmentation permutation reproducibility.
#' @param target_slot Character. Target slot for the resulting segment ratios.
#' @param max_breakpoints Numeric. If segmentation method is CBS performs a pre-
#' segmentation filtering to avoid long processing time due to undo.splits = 'prune'
#' @param n_threads Number of threads used to calculate the distance matrix.
#' Passed to `parallel::mclapply`. As default it uses 1/4 of the detected cores available.
#'
#' @return The segment profile for all cells inside the scCNA object. Can be retrieved with \code{copykit::segment_ratios()}
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @importFrom dplyr mutate bind_cols
#' @importFrom stringr str_detect str_remove str_replace
#' @importFrom S4Vectors metadata
#' @importMethodsFrom SummarizedExperiment assay
#' @export
#'
#' @examples
runSegmentation <- function(scCNA,
                            method = "CBS",
                            genome = "hg38",
                            seed = 17,
                            target_slot = 'segment_ratios',
                            max_breakpoints = 60,
                            n_threads = parallel::detectCores() / 4) {

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Thu Apr  8 16:01:45 2021
  # Checks
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Thu Apr  8 16:01:50 2021

  # check for genome assembly
  if (genome %!in% c("hg19", "hg38")) {
    stop("Genome assembly must be 'hg19' or 'hg38'")
  }

  if (genome != metadata(scCNA)$genome) {
    stop(paste("Incompatible genome assembly, scCNA object was created with",
               metadata(scCNA)$genome, "and runSegmentation is set to:", genome))
  }

  message(paste0("Running segmentation algorithm: ", method, " for genome ", genome))
  message(paste0("Using ", n_threads, " cores."))
  message("Imagine a progress bar here ...")

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

  # genome assembly
  # Reading hg38 VarBin ranges
  if (genome == "hg38") {

    hg38_rg_mod <- hg38_rg
    #match for chrY presence
    chr_sccna <- as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
    hg38_rg_mod <- hg38_rg_mod[which(hg38_rg_mod$chr %in% chr_sccna),]

    hg38_rg_mod <- hg38_rg_mod %>%
      dplyr::mutate(chr = stringr::str_replace(chr, "X", "23"),
                    chr = stringr::str_replace(chr, "Y", "24"))

    chr_info <-  as.numeric(stringr::str_remove(hg38_rg_mod$chr, "chr"))

    ref <- hg38_rg_mod

  }

  # reading hg19 varbin ranges
  if (genome == "hg19") {

    hg19_rg_mod <- hg19_rg
    #match for chrY presence
    chr_sccna <- as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
    hg19_rg_mod <- hg19_rg_mod[which(hg19_rg_mod$chr %in% chr_sccna),]

    hg19_rg_mod <- hg19_rg_mod %>%
      dplyr::mutate(chr = stringr::str_replace(chr, "X", "23"),
                    chr = stringr::str_replace(chr, "Y", "24"))

    chr_info <-  as.numeric(stringr::str_remove(hg19_rg_mod$chr, "chr"))

    ref <- hg19_rg_mod

  }

  if (method == "CBS") {

    if (S4Vectors::metadata(scCNA)$vst == 'ft') {

      counts_df <- assay(scCNA, 'ft')

      CBS_seg <- parallel::mclapply(counts_df, function(x) {
        CNA_object <-
          DNAcopy::CNA(
            x,
            chr_info,
            ref$start,
            data.type = "logratio",
            sampleid = names(x)
          )
        set.seed(seed)
        smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
        segment_smoothed_CNA_object <-
          DNAcopy::segment(
            smoothed_CNA_object,
            alpha = 0.01,
            min.width = 5,
            undo.splits = "prune"
          )
        short_cbs <- segment_smoothed_CNA_object[[2]]
        log_seg_mean_LOWESS <-
          rep(short_cbs$seg.mean, short_cbs$num.mark)
        merge_obj <-
          .MergeLevels(smoothed_CNA_object[, 3], log_seg_mean_LOWESS)$vecMerged
        merge_ratio <- .invft(merge_obj)

      }, mc.cores = n_threads)

    }

    if (S4Vectors::metadata(scCNA)$vst == 'log') {

      # Pre-filtering step for CBS.
      # undo.splits = prune performs better undo of breakpoints without removing
      # true breakpoints. However, in cell with large number of breakpoints it can
      # take more than one hour to finish. The pre-filtering step uses undo.splits = 'sdundo'
      # and removes cells with breakpoints more than the maximum number of breakpoints tolerated
      # this happens in cells that have very noise bin counts profiles.

      message('Performing segmentation pre-filtering.')
      counts_df_pf <- assay(scCNA, 'log')

      CBS_seg_pf <- parallel::mclapply(as.data.frame(counts_df_pf), function(x) {
        CNA_object <-
          DNAcopy::CNA(
            x,
            chr_info,
            ref$start,
            data.type = "logratio",
            sampleid = names(x)
          )
        set.seed(seed)
        smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
        segment_smoothed_CNA_object <-
          DNAcopy::segment(
            smoothed_CNA_object,
            alpha = 0.01,
            min.width = 5,
            undo.splits = "sdundo"
          )
        short_cbs <- segment_smoothed_CNA_object[[2]]
        log_seg_mean_LOWESS <-
          rep(short_cbs$seg.mean, short_cbs$num.mark)
        merge_obj <-
          .MergeLevels(smoothed_CNA_object[, 3], log_seg_mean_LOWESS)$vecMerged
        merge_ratio <- 2^merge_obj

      }, mc.cores = n_threads)

      cbs_seg_df_pf <- dplyr::bind_cols(CBS_seg_pf) %>%
        as.data.frame()

      # copy of scCNA object that later will be removed
      scCNA_pf <- scCNA

      SummarizedExperiment::assay(scCNA_pf, 'segment_ratios') <-
        apply(cbs_seg_df_pf, 2, function(x) x/mean(x)) %>%
        as.data.frame()

      scCNA_pf <- .countBreakpoints(scCNA_pf)

      n_bkpt_cell_remove <- SummarizedExperiment::colData(scCNA_pf) %>%
        as.data.frame() %>%
        dplyr::select(sample, breakpoint_count) %>%
        dplyr::filter(breakpoint_count > max_breakpoints) %>%
        dplyr::pull(sample)

      scCNA <- scCNA[,SummarizedExperiment::colData(scCNA)$sample %!in% n_bkpt_cell_remove]

      rm(scCNA_pf)
      message("Finished pre-filtering.")

      # segmentation with undo.splits = "prune"

      counts_df <- assay(scCNA, 'log')

      CBS_seg <- parallel::mclapply(as.data.frame(counts_df), function(x) {
        CNA_object <-
          DNAcopy::CNA(
            x,
            chr_info,
            ref$start,
            data.type = "logratio",
            sampleid = names(x)
          )
        set.seed(seed)
        smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
        segment_smoothed_CNA_object <-
          DNAcopy::segment(
            smoothed_CNA_object,
            alpha = 0.01,
            min.width = 5,
            undo.splits = "prune"
          )
        short_cbs <- segment_smoothed_CNA_object[[2]]
        log_seg_mean_LOWESS <-
          rep(short_cbs$seg.mean, short_cbs$num.mark)
        merge_obj <-
          .MergeLevels(smoothed_CNA_object[, 3], log_seg_mean_LOWESS)$vecMerged
        merge_ratio <- 2^merge_obj

      }, mc.cores = n_threads)

    }

    cbs_seg_df <- dplyr::bind_cols(CBS_seg) %>%
      as.data.frame()

    # calculating ratios
    scCNA <- calcRatios(scCNA, assay = 'bin_counts')

    SummarizedExperiment::assay(scCNA, target_slot) <-
      apply(cbs_seg_df, 2, function(x) x/mean(x)) %>%
      as.data.frame()

    message("Done.")

    return(scCNA)

  }

  if (method == 'WBS') {

    counts_df <- assay(scCNA, 'log')

    WBS_seg <- parallel::mclapply(as.data.frame(counts_df), function(i) {
      seg <- wbs::wbs(i)
      seg_means <- wbs::means.between.cpt(seg$x,
                                          changepoints(seg,penalty="ssic.penalty")$cpt.ic[["ssic.penalty"]])
      seg_means <- 2^(seg_means)
    }, mc.cores = n_threads)

    wbs_seg_df <- dplyr::bind_cols(WBS_seg) %>%
      as.data.frame()

    scCNA <- calcRatios(scCNA, assay = 'bin_counts')

    SummarizedExperiment::assay(scCNA, target_slot) <-
      apply(wbs_seg_df, 2, function(x) x/mean(x)) %>%
      as.data.frame()

    return(scCNA)

  }

}

