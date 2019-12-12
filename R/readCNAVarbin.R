#' Load data from Copy Number Experiment
#'
#' Creates an S4 class scCNA from the output directory of the copy number pipeline.
#'
#' @param dir A path for the output of the copy number pipeline.
#' @param remove_Y (default == FALSE) If set to TRUE, removes information from the chrY from the dataset.
#'
#' @return A S4 class object containing the segment ratios for each bin (rows) and each sample (columns).
#' @return data can be acessed with \code{segment_ratios} and genomic ranges can be acessed with \code{SummarizedExperiment::rowRanges()}
#'
#' @export
#'
#' @examples
#'
readCNAVarbin <- function(dir,
                          remove_Y = FALSE) {
  # Reads a copy number directory and produces
  # a scCNA object as output

  dat <- suppressMessages(readr::read_tsv(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) %>%
    janitor::clean_names())

  if (remove_Y == TRUE) {
    message("Removed chrY.")
    dat <- dat %>%
      dplyr::filter(chrom != 24)
  }

  #saving data
  seg_data <- dat %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos))

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
    dplyr::filter(row_number() == n()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(end = chrompos + mean_bin_sizes) %>%
    dplyr::pull(end)

  ranges <- rg %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(end = lead(chrompos) - 1) %>%
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
  cna_obj <- scCNA(list(segment_ratios = seg_data),
                   rowRanges = g)

  return(cna_obj)

}
