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
#' @importFrom stringr str_replace
#' @importFrom dplyr rename mutate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
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

  #checks
  if (genome %!in% c("hg19", "hg38")) {
    stop("Genome assembly must be 'hg19' or 'hg38'")
  }

  # Reading hg38 VarBin ranges
  if (genome == "hg38") {

    hg38_rg <- readRDS(here("data/hg38rg.rds"))

    hg38_rg <- hg38_rg %>%
      dplyr::mutate(chr = str_replace(chr, "X", "23"),
             chr = str_replace(chr, "Y", "24"))

    rg <- hg38_rg %>%
      dplyr::rename(Chr = "chr",
                    Start = "start",
                    End = "end") %>%
      dplyr::mutate(GeneID = 1:nrow(hg38_rg))

  }

  # reading hg19 varbin ranges
  if (genome == "hg19") {

    hg19_rg <- readRDS(here("data/hg19rg.rds"))

    hg19_rg <- hg19_rg %>%
      dplyr::mutate(chr = str_replace(chr, "X", "23"),
             chr = str_replace(chr, "Y", "24"))

    rg <- hg19_rg %>%
      dplyr::mutate(GeneID = 1:nrow(hg19_rg))

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

  #LOWESS GC normalization

  varbin_counts_list_gccor <- parallel::mclapply(varbin_counts_list, function(x) {
    gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
    gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
    exp(log(x) - gc_cor_z$y)*median(x)
  },
  mc.cores = n_threads)

  bam_metrics <- dplyr::bind_cols(varbin_reads_list)

  varbin_counts_df <- dplyr::bind_cols(varbin_counts_list_gccor)
  # filtering low read counts
  varbin_counts_df <- varbin_counts_df[which(colSums(varbin_counts_df) != 0)]

  cna_obj <-  scCNA(
    segment_ratios = matrix(nrow = nrow(varbin_counts_df), ncol = ncol(varbin_counts_df)),
    ratios = matrix(nrow = nrow(varbin_counts_df), ncol = ncol(varbin_counts_df)),
    bin_counts = varbin_counts_df,
    rowRanges = GenomicRanges::makeGRangesFromDataFrame(hg38_rg)
  )

}

