#' Plot gene segment ratios
#'
#' geneCopyPlot is a way to visualize the number of copies for different genes across all the cells. It outputs a violin plot with the desired genes.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param genes vector containing the HUGO Symbol for the genes of interest.
#'
#' @return A violin plot with the segment ratios for the genes of interest.
#' @importFrom BiocGenerics subset
#' @importFrom S4Vectors subjectHits
#'
#' @export
#'
#' @examples
#'


geneCopyPlot <- function(scCNA,
                         genes) {

  # theme setup
  my_theme <- list(
    ggplot2::theme(
      axis.title.x = element_text(colour = "gray28", size = 20),
      axis.text.x = element_text(size = 15),
      axis.title.y = element_text(colour = "gray28", size = 20),
      axis.text.y = element_text(size = 15),
      axis.line.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    ),
    xlab(""),
    ylab("Segment ratio")
  )

  # getting ranges from scCNA object
  ranges <- SummarizedExperiment::rowRanges(scCNA)

  # subsetting to only the desired genes
  # hg19_genes genomic range is saved inside sysdata.rda to avoid loading a lot of annotation packages
  hg19_genes_features <- BiocGenerics:::subset(hg19_genes,
                                SYMBOL %in% genes)

  #finding overlaps
  olaps <-
    suppressWarnings(GenomicRanges::findOverlaps(hg19_genes_features,
                 ranges,
                 ignore.strand = TRUE))

  # creating a data frame that will contain the genes and positions (index) in the
  # pipeline ranges.
  # some genes might overlap more than one range (more than one bin), in this case
  # only one will be kept
  df <-
    tibble::tibble(gene = as.character(hg19_genes_features$SYMBOL[S4Vectors::queryHits(olaps)]),
           pos = S4Vectors::subjectHits(olaps),
    ) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  # obtaining seg ratios and sbsetting for the genes
  seg_data <- segment_ratios(scCNA)

  seg_data_genes <- seg_data[df$pos, ] %>%
    dplyr::mutate(gene = df$gene)

  #long format for plotting
  seg_long <- tidyr::gather(seg_data_genes,
                            key = "sample",
                            value = "segratio",
                            -gene)
  #plotting
  p <- ggplot2::ggplot(seg_long, aes(x = gene,
                       y = segratio + 1e-3)) +
    ggplot2::geom_violin() +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    my_theme

  print(p)

}

