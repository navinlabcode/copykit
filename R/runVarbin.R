#' Read BAM file and perform varbin analysis
#'
#' runVarbin performs the variable binning (VarBin) algorithm to a series of BAM files resulting from short-read sequencing.
#'
#' @author Darlan Conterno Minussi
#'
#' @param dir A path for the directory containing .BAM files from short-read sequencing.
#' @param genome Name of the genome assembly. Default: 'hg38'.
#' @param bin_size The resolution of the VarBin method. Default: '200kb'.
#' @param remove_Y (default == FALSE) If set to TRUE, removes information from the chrY from the dataset.
#' @param n_threads Number of cores used to process files
#'
#' @return Segment ratios can be acessed with \code{copykit::segment_ratios}.
#' @return Ratios can be acessed with \code{copykit::ratios}.
#' @return Bin counts can be acessed with \code{copykit::bin_counts}.
#' @return Genomic ranges can be acessed with \code{SummarizedExperiment::rowRanges()}
#'
#' @importFrom Rsubread featureCounts
#' @importFrom stringr str_replace str_remove str_detect
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
                      remove_Y = FALSE,
                      n_threads =  parallel::detectCores() / 4) {


  #checks
  if (genome %!in% c("hg19", "hg38")) {
    stop("Genome assembly must be 'hg19' or 'hg38'")
  }

  # Reading hg38 VarBin ranges
  if (genome == "hg38") {

    rg <- hg38_rg %>%
      dplyr::rename(Chr = "chr",
                    Start = "start",
                    End = "end") %>%
      dplyr::mutate(GeneID = 1:nrow(hg38_rg))

    if (remove_Y == TRUE) {

      rg <- dplyr::filter(rg,
                          Chr != "chrY")

    }

  }

  # reading hg19 varbin ranges
  if (genome == "hg19") {

    rg <- hg19_rg %>%
      dplyr::mutate(GeneID = 1:nrow(hg19_rg))

    if (remove_Y == TRUE) {

      rg <- dplyr::filter(rg,
                          Chr != "chrY")

    }

  }

  files <- list.files(dir, pattern = "*.bam", full.names = T)

  # managing .bai files
  if (any(sapply(files, function(x) !stringr::str_detect(x, ".bai")))) {
    files <- files[!stringr::str_detect(files, ".bai")]
  }

  if (any(sapply(files, function(x) !stringr::str_detect(x, ".bam")))) {
    stop("Directory does not contain .bam files.")
  }

  files_names <- list.files(dir, pattern = "*.bam", full.names = F)
  files_names <- stringr::str_remove(files_names, ".bam")

  # managing .bai files on filenames
  if (any(sapply(files_names, function(x) !stringr::str_detect(x, ".bai")))) {
    files_names <- files_names[!stringr::str_detect(files_names, ".bai")]
  }

  varbin_counts_list_all_fields <- parallel::mclapply(files, Rsubread::featureCounts,
                                                      ignoreDup = TRUE,
                                                      countMultiMappingReads = FALSE,
                                                      annot.ext = rg,
                                                      useMetaFeatures = FALSE,
                                                      mc.cores = n_threads)

  varbin_counts_list <- lapply(varbin_counts_list_all_fields,
                               '[[',
                               1)



  names(varbin_counts_list) <- files_names

browser()
  #LOWESS GC normalization

  varbin_counts_list_gccor <- parallel::mclapply(varbin_counts_list, function(x) {
    gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
    gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
    exp(log(x) - gc_cor_z$y)*median(x)
  },
  mc.cores = n_threads)

  varbin_counts_df <- dplyr::bind_cols(varbin_counts_list_gccor)
  # filtering low read counts
  varbin_counts_df <- varbin_counts_df[which(colSums(varbin_counts_df) != 0)]

  #sample name to metadata
  SummarizedExperiment::colData(cna_obj)$sample <- names(varbin_counts_df)
  colnames(cna_obj) <- names(varbin_counts_df)


  rg <- rg %>%
    dplyr::select(-strand,
                  -GeneID)

  rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,
                                                   ignore.strand = TRUE,
                                                   keep.extra.columns = TRUE)


  cna_obj <-  scCNA(
    segment_ratios = as.data.frame(matrix(nrow = nrow(varbin_counts_df), ncol = ncol(varbin_counts_df))),
    ratios = as.data.frame(matrix(nrow = nrow(varbin_counts_df), ncol = ncol(varbin_counts_df))),
    bin_counts = as.data.frame(varbin_counts_df),
    rowRanges = rg_gr
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
  # ADDING READS METRICS TO METADATA
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021

  varbin_reads_list <- lapply(varbin_counts_list_all_fields,
                              '[[',
                              4)

  names(varbin_reads_list) <- files_names

  # saving info and removing columns from list elements
  metadata_info_names <- varbin_reads_list[[1]][,1]
  varbin_reads_info <- lapply(varbin_reads_list, function(x)  x[,-1, drop = FALSE])

  bam_metrics <- cbind(metadata_info_names, dplyr::bind_cols(varbin_reads_info))



  return(cna_obj)


}

