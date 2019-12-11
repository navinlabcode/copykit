#' Load data from Copy Number Experiment
#'
#' Creates an S4 class scCNA from the output directory of the copy number pipeline.
#'
#' @param dir A path for the output of the copy number pipeline.
#'
#' @return A S4 class object containing the segment ratios for each bin (rows) and each sample (columns).
#' @export
#'
#' @examples
readCNAVarbin <- function(dir) {
  # Reads a copy number directory and produces
  # a scCNA object as output

  seg_data <- readr::read_tsv(fs::dir_ls(
    path = dir,
    recurse = T,
    glob = "*uber*seg.txt"
  )) %>%
    janitor::clean_names() %>%
    dplyr::select(-c(chrom,
                     chrompos,
                     abspos)) %>%
    as.data.frame()

  cna_obj <- scCNA(list(segment_ratios = seg_data))

  return(cna_obj)

}
