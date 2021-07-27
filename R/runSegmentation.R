#' Run Segmentation
#'
#' Runs a segmentation algorithm using the ratio data.
#'
#' @param scCNA The scCNA object
#' @param method Character. Segmentation method of choice.
#' @param seed Numeric. Set seed for CBS segmentation permutation reproducibility
#' @param undo.splits A character string specifying how change-points are to be
#' undone, if at all. Default is "none". Other choices are "prune", which uses
#'  a sum of squares criterion, and "sdundo", which undoes splits that are not
#'  at least this many SDs apart. See \code{\link[DNAcopy]{segment}}
#' @param name Character. Target slot for the resulting segment ratios.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @return The segment profile for all cells inside the scCNA object.
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @importFrom dplyr mutate bind_cols
#' @importFrom stringr str_detect str_remove str_replace
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay
#' @importFrom BiocParallel bplapply bpparam
#' @importMethodsFrom SummarizedExperiment assay
#' @export
#'
#' @examples
runSegmentation <- function(scCNA,
                            method = "CBS",
                            seed = 17,
                            undo.splits = 'sdundo',
                            name = 'segment_ratios',
                            BPPARAM = bpparam()) {
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Thu Apr  8 16:01:45 2021
  # Checks
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Thu Apr  8 16:01:50 2021
  # genome assembly
  if (S4Vectors::metadata(scCNA)$genome == "hg19") {
    genome <- "hg19"
  }

  if (S4Vectors::metadata(scCNA)$genome == "hg38") {
    genome <- "hg38"
  }

  message(paste0(
    "Running segmentation algorithm: ",
    method,
    " for genome ",
    genome
  ))

  # genome assembly
  # Reading hg38 VarBin ranges
  if (genome == "hg38") {
    hg38_rg_mod <- hg38_rg
    #match for chrY presence
    chr_sccna <-
      as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
    hg38_rg_mod <-
      hg38_rg_mod[which(hg38_rg_mod$chr %in% chr_sccna), ]

    hg38_rg_mod <- hg38_rg_mod %>%
      dplyr::mutate(
        chr = stringr::str_replace(chr, "X", "23"),
        chr = stringr::str_replace(chr, "Y", "24")
      )

    chr_info <-
      as.numeric(stringr::str_remove(hg38_rg_mod$chr, "chr"))

    ref <- hg38_rg_mod

  }

  # reading hg19 varbin ranges
  if (genome == "hg19") {
    hg19_rg_mod <- hg19_rg
    #match for chrY presence
    chr_sccna <-
      as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
    hg19_rg_mod <-
      hg19_rg_mod[which(hg19_rg_mod$chr %in% chr_sccna), ]

    hg19_rg_mod <- hg19_rg_mod %>%
      dplyr::mutate(
        chr = stringr::str_replace(chr, "X", "23"),
        chr = stringr::str_replace(chr, "Y", "24")
      )

    chr_info <-
      as.numeric(stringr::str_remove(hg19_rg_mod$chr, "chr"))

    ref <- hg19_rg_mod

  }



  if (method == "CBS") {

    if (S4Vectors::metadata(scCNA)$vst == 'ft') {
      counts_df <- SummarizedExperiment::assay(scCNA, 'ft')

      CBS_seg <- BiocParallel::bplapply(counts_df, FUN = function(x) {
        CNA_object <-
          DNAcopy::CNA(x,
                       chr_info,
                       ref$start,
                       data.type = "logratio",
                       sampleid = names(x))
        set.seed(seed)
        smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
        segment_smoothed_CNA_object <-
          .quiet(
            DNAcopy::segment(
              smoothed_CNA_object,
              alpha = 0.01,
              min.width = 5,
              undo.splits = undo.splits
            )
          )
        short_cbs <- segment_smoothed_CNA_object[[2]]
        log_seg_mean_LOWESS <-
          rep(short_cbs$seg.mean, short_cbs$num.mark)
        merge_obj <-
          .MergeLevels(smoothed_CNA_object[, 3], log_seg_mean_LOWESS)$vecMerged
        merge_ratio <- .invft(merge_obj)

      }, BPPARAM = BPPARAM)

    }

    if (S4Vectors::metadata(scCNA)$vst == 'log') {

      warning("undo.splits = 'prune' and 'log' assay can run for long time
              for cells with large number of breakpoints",
              noBreaks. = TRUE)
      # segmentation with undo.splits = "prune"

      counts_df <- SummarizedExperiment::assay(scCNA, 'log')

      CBS_seg <-
        BiocParallel::bplapply(as.data.frame(counts_df), function(x) {
          CNA_object <-
            DNAcopy::CNA(x,
                         chr_info,
                         ref$start,
                         data.type = "logratio",
                         sampleid = names(x))
          set.seed(seed)
          smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
          segment_smoothed_CNA_object <-
            .quiet(
              DNAcopy::segment(
                smoothed_CNA_object,
                alpha = 0.01,
                min.width = 5,
                undo.splits = undo.splits
              )
            )
          short_cbs <- segment_smoothed_CNA_object[[2]]
          log_seg_mean_LOWESS <-
            rep(short_cbs$seg.mean, short_cbs$num.mark)
          merge_obj <-
            .MergeLevels(smoothed_CNA_object[, 3], log_seg_mean_LOWESS)$vecMerged
          merge_ratio <- 2 ^ merge_obj


        }, BPPARAM = BPPARAM)

    }

    cbs_seg_df <- dplyr::bind_cols(CBS_seg) %>%
      as.data.frame()

    # calculating ratios
    scCNA <- calcRatios(scCNA, assay = 'bin_counts')

    SummarizedExperiment::assay(scCNA, name) <-
      apply(cbs_seg_df, 2, function(x)
        x / mean(x)) %>%
      as.data.frame()

    message("Done.")

    return(scCNA)

  }

}
