#' Load data from Copy Number Experiment
#'
#' Creates an S4 class scCNA from the output directory of the copy number pipeline.
#' readVarbinCNA searches for the uber*.seg uber.bin and uber.ratio files resulting from the Varbin copy number pipeline in the provided directory.
#' The scCNA object contains the segment ratios, ratios and bincounts within the assay slot. where each bin is row and each sample (cell) is a column.
#' Genomic ranges are stored in a GRanges object containing chromosome number, start coordinate, end cordinate and absolute genomic position. Each row represents the coordinates for one bin.
#'
#' @author Darlan Conterno Minussi
#'
#' @param dir A path for the output of the copy number pipeline.
#' @param remove_Y (default == FALSE) If set to TRUE, removes information from the chrY from the dataset.
#'
#' @return Segment ratios can be acessed with \code{copykit::segment_ratios}.
#' @return Ratios can be acessed with \code{copykit::ratios}.
#' @return Bin counts can be acessed with \code{copykit::bin_counts}.
#' @return Genomic ranges can be acessed with \code{SummarizedExperiment::rowRanges()}
#'
#' @export
#'
#' @examples
#'
readVarbinCNA <- function(dir,
                          remove_Y = FALSE) {
  # Reads a copy number directory and produces
  # a scCNA object as output

  # checks
  if (fs::file_exists(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) == FALSE) {
    stop(
      "Segment ratio matrix can't be found in the provided directory. Please make sure uber.seg file can be found."
    )
  }

  if (length(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) > 1) {
    stop(
      "More than one uber.seg file can be found at the provided directory. Please make sure to only have one sample at that location."
    )
  }

  # reading segment ratios
  message("Importing segment ratios.")
  dat <- readr::read_tsv(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  ),
  col_types = readr::cols()) %>%
    janitor::clean_names()

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
  dat_rat <- readr::read_tsv(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*ratio.txt"
  ),
  col_types = readr::cols()) %>%
    janitor::clean_names()

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
  dat_bin <- readr::read_tsv(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*bin.txt"
  ),
  col_types = readr::cols()) %>%
    janitor::clean_names()

  if (remove_Y == TRUE) {
    dat_bin <- dat_bin %>%
      dplyr::filter(chrom != 24)
  }

  dat_bin <- dat_bin %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))

  if (remove_Y == TRUE) {
    message("Removed ChrY information.")
  }

  # saving ranges, generating start and end
  # Unfortunately, varbin does not output start and end position, so it has to be derived from the data
  # leaving a problem to the correct end value of the bin. We can't use the chromosome length due to the
  # telomeric regions being blacklisted so I opted to add the average of the bin size to the last bin
  # this will be solved by adding start and end to the varbin pipeline but won't be backwards compatible with previously
  # processed datasets

  rg <- dat %>%
    dplyr::select(c(chrom,
                    chrompos,
                    abspos)) %>%
    as.data.frame()

  mean_bin_sizes <- rg %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(bin_size = (dplyr::lead(chrompos) - 1) - chrompos) %>%
    tidyr::drop_na() %>%
    dplyr::summarise(mean_bin_size = mean(bin_size)) %>%
    dplyr::pull(mean_bin_size)

  ranges_end <- rg %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(end = dplyr::lead(chrompos) - 1) %>%
    dplyr::filter(dplyr::row_number() == dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(end = chrompos + mean_bin_sizes) %>%
    dplyr::pull(end)

  ranges <- rg %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(end = dplyr::lead(chrompos) - 1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      chr = chrom,
      start = chrompos,
      end = end,
      abspos = abspos
    ) %>%
    dplyr::select(chr,
                  start,
                  end,
                  abspos)

  if (length(which(is.na(ranges$end) == TRUE)) == length(ranges_end)) {
    ranges$end[is.na(ranges$end)] <- ranges_end
  } else
    stop("Problem with adding the value of the last bin end")

  g <- GenomicRanges::makeGRangesFromDataFrame(ranges,
                                               keep.extra.columns = TRUE,
                                               ignore.strand = T)


  # creating scCNA object
  cna_obj <- scCNA(
    segment_ratios = seg_data,
    ratios = dat_rat,
    bin_counts = dat_bin,
    rowRanges = g
  )

  #sample name to metadata
  SummarizedExperiment::colData(cna_obj)$sample <- names(seg_data)

  return(cna_obj)

}
