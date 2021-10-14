#' Run Segmentation
#'
#' Runs a segmentation algorithm using the ratio data.
#'
#' @param scCNA The scCNA object
#' @param assay String with the name of the assay to pull data from to run the
#' segmentation.
#' @param method A character with the segmentation method of choice.
#' @param alpha A numeric with the. significance levels for the test to accept
#' change-points for CBS segmentation. See \code{\link[DNAcopy]{segment}}.
#' @param gamma A numeric passed on to 'multipcf' segmentation. Penalty for each
#'  discontinuity in the curve, default is 40. See \code{\link[copynumber]{multipcf}}.
#' @param seed Numeric. Set seed for CBS segmentation permutation reproducibility
#' @param undo.splits A character string specifying how change-points are to be
#' undone, if at all. Default is "none". Other choices are "prune", which uses
#'  a sum of squares criterion, and "sdundo", which undoes splits that are not
#'  at least this many SDs apart. See \code{\link[DNAcopy]{segment}}
#' @param name Character. Target slot for the resulting segment ratios.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details
#'
#' \itemize{
#'
#'    \item{CBS:} #' \code{runSegmentation} Fits a piece-wise constant function
#'    to the transformed the smoothed bin counts. Bin counts are smoothed with
#'    \code{\link[DNAcopy]{smooth.CNA}} using the Circular Binary Segmentation
#'    (CBS) algorithm from \code{\link[DNAcopy]{segment}} with default it applies
#'    undo.prune with value of 0.05.
#'
#'    \item{multipcf:} Performs the joint segmentation from the \code{copynumber}
#'    package using the \code{\link[copynumber]{multipcf}} function. By fitting
#'    piecewise constant curves with common breakpoints for all samples.
#'
#'
#' }
#'
#'
#' The resulting segment means are further refined with MergeLevels to join
#' adjacent segments with non-significant differences in segmented means.
#'
#' @return The segment profile for all cells inside the scCNA object.
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @importFrom aCGH mergeLevels combine.func
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
                            method = c("CBS", "multipcf"),
                            seed = 17,
                            alpha = 1e-9,
                            gamma = 40,
                            undo.splits = 'prune',
                            name = 'segment_ratios',
                            BPPARAM = bpparam()) {
  # Args
  method <- match.arg(method)

  #checks
  if (!is.numeric(alpha)) stop("Argument alpha must be numeric.")
  if (!is.numeric(gamma)) stop("Argument gamma must be numeric.")
  if (!is.numeric(seed)) stop("Argument seed must be numeric.")

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Genome Assembly
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # adding genome assembly info to metadata
  if (S4Vectors::metadata(scCNA)$genome == "hg19") {
    genome <- "hg19"
  }

  if (S4Vectors::metadata(scCNA)$genome == "hg38") {
    genome <- "hg38"
  }

  # Reading hg38 VarBin ranges
  if (genome == "hg38") {
    hg38_rg_mod <- hg38_rg
    #match for chrY presence
    chr_sccna <-
      as.character(as.data.frame(SummarizedExperiment::rowRanges(scCNA))$seqnames)
    hg38_rg_mod <-
      hg38_rg_mod[which(hg38_rg_mod$chr %in% chr_sccna),]

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
      hg19_rg_mod[which(hg19_rg_mod$chr %in% chr_sccna),]

    hg19_rg_mod <- hg19_rg_mod %>%
      dplyr::mutate(
        chr = stringr::str_replace(chr, "X", "23"),
        chr = stringr::str_replace(chr, "Y", "24")
      )

    chr_info <-
      as.numeric(stringr::str_remove(hg19_rg_mod$chr, "chr"))

    ref <- hg19_rg_mod

  }

  ref_chrarm <- ref %>%
    dplyr::mutate(chrarm = paste0(stringr::str_remove(chr, 'chr'), arm)) %>%
    dplyr::mutate(chrarm = chrarm)

  levels_chrarm <- gtools::mixedsort(unique(ref_chrarm$chrarm))

  ref_chrarm <- ref_chrarm %>%
    dplyr::mutate(chrarm = as.factor(chrarm)) %>%
    dplyr::mutate(chrarm = forcats::fct_relevel(chrarm, levels_chrarm))

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Data Setup
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (S4Vectors::metadata(scCNA)$vst == 'ft') {
    counts_df <- SummarizedExperiment::assay(scCNA, 'ft')
  }

  if (S4Vectors::metadata(scCNA)$vst == 'log') {
    counts_df <- SummarizedExperiment::assay(scCNA, 'log')
  }

  # smoothing data
  message("Smoothing bin counts.")
  smooth_counts <-
    BiocParallel::bplapply(as.data.frame(counts_df), function(x) {
      CNA_object <-
        DNAcopy::CNA(x,
                     ref_chrarm$chrarm,
                     ref$start,
                     data.type = "logratio",
                     sampleid = names(x))
      set.seed(seed)
      smoothed_CNA_counts <- DNAcopy::smooth.CNA(CNA_object)[, 3]
    }, BPPARAM = BPPARAM)

  smooth_counts_df <- dplyr::bind_cols(smooth_counts) %>%
    as.data.frame() %>%
    round(2)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Segmentation methods
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  message(paste0(
    "Running segmentation algorithm: ",
    method,
    " for genome ",
    genome
  ))

  if (method == "CBS") {
    seg_list <-
      BiocParallel::bplapply(
        smooth_counts_df,
        FUN = function(x) {
          CNA_object <-
            DNAcopy::CNA(x,
                         ref_chrarm$chrarm,
                         ref$start,
                         data.type = "logratio",
                         sampleid = names(x))
          set.seed(seed)

          segment_smoothed_CNA_object <-
            .quiet(
              DNAcopy::segment(
                CNA_object,
                alpha = alpha,
                min.width = 5,
                undo.splits = undo.splits
              )
            )
          short_cbs <- segment_smoothed_CNA_object[[2]]
          log_seg_mean_LOWESS <-
            rep(short_cbs$seg.mean, short_cbs$num.mark)

        },
        BPPARAM = BPPARAM
      )

    seg_df <- dplyr::bind_cols(seg_list) %>%
      as.data.frame() %>%
      round(2)

  }



  if (method == "multipcf") {

    smooth_multipcf <- cbind(as.numeric(str_remove(ref_chrarm$chr, 'chr')),
                             ref_chrarm$start,
                             smooth_counts_df)

    mpcf <- copynumber::multipcf(smooth_multipcf,
                         arms = str_extract(ref_chrarm$chrarm,"[pq]"))

    seg_df <- apply(mpcf[,6:ncol(mpcf)], 2, function(x) {
      rep.int(x, mpcf$n.probes)
    })

    seg_df <- round(as.data.frame(seg_df),2)


  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # merge levels
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  message("Merging levels.")

  if (S4Vectors::metadata(scCNA)$vst == 'ft') {
    smooth_counts_df[smooth_counts_df == 0] <- 1e-4
    seg_df[seg_df == 0] <- 1e-4
    smooth_counts_df <- log(smooth_counts_df)
    seg_df <- log(seg_df)
  }

  seg_ml_list <- BiocParallel::bplapply(seq_along(seg_df), function(i) {

    cell_name <- names(seg_df)[i]
    smoothed_cell_ct <- smooth_counts_df[,i]
    seg_means_cell <- seg_df[,i]
    seg_means_ml <- aCGH::mergeLevels(smoothed_cell_ct,
                                      seg_means_cell,
                                      verbose = 0,
                                      pv.thres = 1e-10)$vecMerged
  })

  names(seg_ml_list) <- names(seg_df)
  seg_ml_df <- dplyr::bind_cols(seg_ml_list)

  if (S4Vectors::metadata(scCNA)$vst == 'ft') {
    seg_ml_df <- round(2^seg_ml_df,2)
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # reverting the transformation back to ratios
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (S4Vectors::metadata(scCNA)$vst == 'ft') {

    SummarizedExperiment::assay(scCNA, 'smoothed_bincounts') <-
      .invft(2^smooth_counts_df)

    seg_ratio_df <- .invft(seg_ml_df)
  }

  if (S4Vectors::metadata(scCNA)$vst == 'log') {

    SummarizedExperiment::assay(scCNA, 'smoothed_bincounts') <-
      2^smooth_counts_df

    seg_ratio_df <- 2 ^ seg_ml_df
  }

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # saving information to the scCNA object
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #saving as segment ratios
  seg_ratios <- sweep(seg_ratio_df, 2, apply(seg_ratio_df, 2, mean), '/')
  SummarizedExperiment::assay(scCNA, name) <- round(seg_ratios,2)


  # calculating ratios from the bincounts, used for ratio plots
  scCNA <- calcRatios(scCNA, assay = 'smoothed_bincounts')

  message("Done.")

  return(scCNA)

}
