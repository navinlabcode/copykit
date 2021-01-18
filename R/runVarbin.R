#' Read BAM file and perform varbin analysis
#'
#' runVarbin performs the variable binning (VarBin) algorithm to a series of BAM files resulting from short-read sequencing.
#'
#' @author Darlan Conterno Minussi
#'
#' @param dir A path for the directory containing .BAM files from short-read sequencing.
#' @param genome Name of the genome assembly. Default: 'hg38'.
#' @param bin_size The resolution of the VarBin method. Default: '200kb'.
#' @param n_threads Number of cores used to process files
#'
#' @return Segment ratios can be acessed with \code{copykit::segment_ratios}.
#' @return Ratios can be acessed with \code{copykit::ratios}.
#' @return Bin counts can be acessed with \code{copykit::bin_counts}.
#' @return Genomic ranges can be acessed with \code{SummarizedExperiment::rowRanges()}
#'
#' @importFrom Rsubread featureCounts
#'
#' @export
#'
#' @examples
#'
#'

runVarbin <- function(dir,
                      genome = "hg38",
                      bin_size = "200kb",
                      n_threads = 1 ) {

  # Reading hg38 VarBin ranges
  if (geneome == "hg38") {

    hg38_rg <- readRDS(here("data/hg38rg.rds"))

    hg38_rg <- hg38_rg %>%
      mutate(chr = str_replace(chr, "X", "23"),
             chr = str_replace(chr, "Y", "24"))

    rg <- hg38_rg %>%
      dplyr::rename(Chr = "chr",
                    Start = "start",
                    End = "end") %>%
      mutate(GeneID = 1:nrow(hg38_rg))

  }
   dir = "/volumes/seq/projects/CNA_projects/DT_CNA/snap_frozen/Breast/TNBC/TN20/TN20_2020_06_19_combine_new_run_output/output/sort/"

  files <- list.files(dir, pattern = "*.bam", full.names = T)[1:10]

  if (!str_detect(files, ".bam")) {
    stop("Directory does not contain .bam files.")
  }

  files_names <- list.files(dir, pattern = "*.bam", full.names = F)[1:10]
  files_names <- stringr::str_remove(files_names, ".bam")

  varbin_counts_list_all_fields <- parallel::mclapply(files, Rsubread::featureCounts,
                                                      ignoreDup = TRUE,
                                                      countMultiMappingReads = FALSE,
                                                      annot.ext = rg,
                                                      useMetaFeatures = FALSE,
                                                      mc.cores = n_threads)

  varbin_counts_list <- lapply(varbin_counts_list_all_fields,
                               '[[',
                               1)

  varbin_reads_list <- lapply(varbin_counts_list_all_fields,
                               '[[',
                               4)

  names(varbin_counts_list) <- files_names
  names(varbin_reads_list) <- files_names

  varbin_reads_list <- lapply(varbin_reads_list)

  # filtering low read counts
  bam_metrics <- dplyr::bind_cols(varbin_reads_list)


  varbin_counts_df <- dplyr::bind_cols(varbin_counts_list)

  # segmenting
  ratios_df <- sweep(varbin_counts_df, 2, apply(varbin_counts_df, 2, median), '/')

  CBS_seg <- parallel::mclapply(ratios_df, function(x) {
    CNA_object <- DNAcopy::CNA(log(x+1e-3, base=2), as.numeric(stringr::str_remove(hg38_rg$chr, "chr")), hg38_rg$start, data.type="logratio",
                 sampleid=names(x))
    smoothed_CNA_object <- DNAcopy::smooth.CNA(CNA_object)
    segment_smoothed_CNA_object <- DNAcopy::segment(smoothed_CNA_object, alpha=0.01, min.width=5, undo.splits="prune", undo.prune=0.05)
    short_cbs <- segment_smoothed_CNA_object[[2]]
    log_seg_mean_LOWESS <- rep(short_cbs$seg.mean, short_cbs$num.mark)
    merge_obj <- MergeLevels(smoothed_CNA_object[,3],log_seg_mean_LOWESS)$vecMerged

  }, mc.cores = n_threads
)

  cbs_seg_df <- bind_cols(CBS_seg) %>%
    as.data.frame()

  cna_obj <-   cna_obj <- scCNA(
    segment_ratios = seg_data,
    ratios = dat_rat,
    bin_counts = varbin_counts_df,
    rowRanges = g
  )

}

