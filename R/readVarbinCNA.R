#' Load data from Copy Number Experiment
#'
#' Creates an S4 class scCNA from the output directory of the copy number pipeline.
#' readVarbinCNA searches for the uber*.seg uber.bin and uber.ratio files resulting from the Varbin copy number pipeline in the provided directory.
#' The scCNA object contains the segment ratios, ratios and bincounts within the assay slot. where each bin is row and each sample (cell) is a column.
#' Genomic ranges are stored in a GRanges object containing chromosome number, start coordinate, end cordinate and absolute genomic position. Each row represents the coordinates for one bin.
#'
#' @author Darlan Conterno Minussi, Yun Yan
#'
#' @param dir A path for the output of the copy number pipeline.
#' @param remove_Y (default == FALSE) If set to TRUE, removes information from the chrY from the dataset.
#' @param genome_version Name of the genome assembly. Default: 'hg19'.
#' @param bin_size The resolution of the VarBin method. Default: '200k'. Available options: '100k', '200k'.
#' @param clean_names Clean sample names using \code{\link{janitor::clean_names}}.
#' @return Segment ratios can be accessed with \code{copykit::segment_ratios}.
#' @return Ratios can be accessed with \code{copykit::ratios}.
#' @return Bin counts can be accessed with \code{copykit::bin_counts}.
#' @return Genomic ranges can be accessed with \code{SummarizedExperiment::rowRanges()}
#'
#' @importFrom dplyr filter select mutate rename
#' @importFrom GenomicRanges seqnames
#' @importFrom data.table fread
#' @importFrom fs dir_ls
#' @import rlang
#' @export
#'
#' @examples
#'
readVarbinCNA <- function(dir,
                          remove_Y = FALSE,
                          genome_version = c('hg19'),
                          bin_size = c('200k', '100k'),
                          clean_names = TRUE) {
  # Reads a copy number directory and produces
  # a scCNA object as output

  # checks
  if (fs::file_exists(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) == FALSE) {
    stop(
      "Segment ratio matrix can't be found in the provided directory.
      Please make sure a uber.seg file can be found."
    )
  }

  # checking for the existence of more than one uber file
  if (length(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) > 1) {
    stop(
      "More than one uber.seg file can be found at the provided directory.
      Please make sure to only have one sample at that location."
    )
  }

  # importing data
  message("Importing segment ratios.")
  dat <- data.table::fread(input = fs::dir_ls(
      path = dir,
      recurse = T,
      glob = "*uber*seg.txt"
    ),
    showProgress = TRUE,
    integer64 = "double") %>%
    as.data.frame()

  if (clean_names == TRUE) {
    dat <- janitor::clean_names(dat)
  }

  if (remove_Y == TRUE) {
    dat <- dat %>%
      dplyr::filter(chrom != 24)
  }

  #saving segment data
  seg_data <- dat %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))


  # reading ratios
  message("Importing ratios.")
  dat_rat <- data.table::fread(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*ratio.txt"
  ),
  showProgress = TRUE,
  integer64 = "double") %>%
    janitor::clean_names() %>%
    as.data.frame()

  if (remove_Y == TRUE) {
    dat_rat <- dat_rat %>%
      dplyr::filter(chrom != 24)
  }

  dat_rat <- dat_rat %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))

  # reading bin counts
  message("Importing bin counts.")
  dat_bin <- data.table::fread(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*bin.txt"
  ),
  showProgress = TRUE,
  integer64 = "double") %>%
    as.data.frame()

  if (clean_names == TRUE) {
    dat_bin <- janitor::clean_names(dat_bin)
  }

  if (remove_Y == TRUE) {
    dat_bin <- dat_bin %>%
      dplyr::filter(chrom != 24)
  }

  dat_bin <- dat_bin %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))

  # Fetch the locations (and other informations) of varbins
  rg <- dat %>%
    dplyr::select(c(chrom,
                    chrompos,
                    abspos)) %>%
    as.data.frame()

  genome_version <- match.arg(genome_version)
  bin_size <- match.arg(bin_size)
  grlist_varbin <- switch(genome_version,
                          hg19 = varbin_hg19_grangeslist)
  tmp_key <- paste0('res_', bin_size)
  gr_varbin_full <- grlist_varbin[[tmp_key]]

  GenomeInfoDb::seqlevelsStyle(gr_varbin_full) <- 'Ensembl'
  gr_varbin_full <- GenomeInfoDb::renameSeqlevels(
    gr_varbin_full, c(X=23, Y=24))

  rg$chrompos <- rg$chrompos + 1  ## To 1-based
  rg$abspos <- rg$abspos + 1  ## To 1-based
  key_query <- paste0(rg$chrom,
                      '_',
                      rg$chrompos)
  key_ref   <- paste0(GenomicRanges::seqnames(gr_varbin_full),
                      '_',
                      IRanges::start(gr_varbin_full))
  idx <- match(key_query, key_ref)
  if (anyNA(idx)) {
    warning('Input ', sum(is.na(idx)),
            'varbins are not recorded in white sheet\n.')
  }
  g <- gr_varbin_full[idx, ]
  g$abspos <- rg$abspos


  g <- GenomeInfoDb::renameSeqlevels(
    g, c(`23`='X', `24`='Y'))
  GenomeInfoDb::seqlevelsStyle(g) <- 'UCSC'  ## add chr prefix
  # A chr1, chr2, ..., chrX, chrY style

  # creating scCNA object
  cna_obj <- scCNA(
    segment_ratios = seg_data,
    ratios = dat_rat,
    bin_counts = dat_bin,
    rowRanges = g
  )

  #sample name to metadata
  SummarizedExperiment::colData(cna_obj)$sample <- names(seg_data)

  # reading metrics
  if (rlang::is_empty(
    fs::dir_ls(path = dir, recurse = T, glob = "*stat_metrics.txt"))){
    warning("No metrics file found. \n
            Metrics are needed if you'd like to run copykit::runMetrics()\n
            Make sure folder metrics with file all_stat_metrics.txt can be found by copykit::runVarbinCNA()")
  } else {
    if (fs::file_exists(
      fs::dir_ls(path = dir, recurse = T, glob = "*stat_metrics.txt"))) {
      message("Importing metrics.")
      dat_metrics <- data.table::fread(
        fs::dir_ls(
          path = dir,
          recurse = T,
          glob = "*stat_metrics.txt"
        ),
        showProgress = TRUE,
        integer64 = "double"
      ) %>%
        janitor::clean_names() %>%
        dplyr::rename(sample = "sample_name") %>%
        as.data.frame()

      if (clean_names == TRUE) {
        dat_metrics <- dat_metrics %>%
          dplyr::mutate(sample = janitor::make_clean_names(sample))
      }

      # adding metrics to metadata
      # making sure they have the same order
      dat_metrics <- dat_metrics[match(SummarizedExperiment::colData(cna_obj)$sample,
                                       dat_metrics$sample),]

      if (identical(dat_metrics$sample,
                    SummarizedExperiment::colData(cna_obj)$sample)) {

        SummarizedExperiment::colData(cna_obj)$reads_total <- dat_metrics$total_reads
        SummarizedExperiment::colData(cna_obj)$reads_assigned_bins <- dat_metrics$reads_kept
        SummarizedExperiment::colData(cna_obj)$percentage_duplicates <- round(dat_metrics$dups_removed/dat_metrics$total_reads,2)

      }

    }
  }

  if (remove_Y == TRUE) {
    message("Removed ChrY information.")
  }

  return(cna_obj)

}
