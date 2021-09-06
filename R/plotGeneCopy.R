#' plotGeneCopy
#'
#' Visualization for gene-wise copy number states
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param genes A vector of strings containing the HUGO Symbol for the gene
#' of interest.
#' @param genome A string with the chosen genome assembly.
#' @param geom A string with the geom for plotting.
#' @param label A string with the name of the column from
#' \code{\link[SummarizedExperiment]{colData}} to color the points
#'
#' @details plotGeneCopy finds overlaps of the varbin scaffolds genomic ranges
#' which can be accessed with \code{\link[SummarizedExperiment]{rowRanges}}
#' with the genes from the assemblies of either hg19 or hg38. The genomic ranges
#' from hg19 comes from package \code{TxDb.Hsapiens.UCSC.hg19.knownGene} whereas
#' for hg38 package \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#'
#' If the argument geom is set to 'barplot' plotGeneCopy calculates the gene-wise
#'  frequencies of each copy number state for the selected genes across all of
#'  the cells. For this reason, geom 'barplot' can only be used with the argument
#'  assay set to 'integer'.
#'
#' @return A ggplot object with a plot of the gene-wise copy number states.
#'
#' @importFrom BiocGenerics subset
#' @importFrom BiocGenerics `%in%`
#' @importFrom S4Vectors subjectHits metadata
#' @importFrom forcats fct_reorder
#' @importFrom SummarizedExperiment assay
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'

plotGeneCopy <- function(scCNA,
                         genes,
                         geom = c("swarm", "barplot", "violin"),
                         label = NULL,
                         assay = "segment_ratios") {

  geom <- match.arg(geom)

  # checks
  # check if label exists
  if (!is.null(label)) {
    metadata <- as.data.frame(SummarizedExperiment::colData(scCNA))

    message(paste0("Coloring by: ", label))
  }

  if (!is.null(label) && !(label %in% colnames(metadata))) {
    stop(paste0("Label ", label, " is not a column of the scCNA object."))
  }

  if (geom == 'barplot' && assay != 'integer') {
    stop("Argument geom 'barplot' can only be used with assay == 'integer'")
  }

  # genome assembly
  if (S4Vectors::metadata(scCNA)$genome == "hg19") {
    genes_assembly <- hg19_genes
  }

  if (S4Vectors::metadata(scCNA)$genome == "hg38") {
    genes_assembly <- hg38_genes
  }

  # theme setup
  my_theme <- list(
    ggplot2::theme(
      axis.title.x = element_text(colour = "gray28", size = 20),
      axis.text.x = element_text(
        size = 15,
        vjust = 0.5,
        hjust = 1,
        angle = 90
      ),
      axis.title.y = element_text(colour = "gray28", size = 20),
      axis.text.y = element_text(size = 15),
      axis.line.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    ),
    xlab(""),
    ylab("segment ratio")
  )

  # getting ranges from scCNA object
  ranges <- SummarizedExperiment::rowRanges(scCNA)

  # subsetting to only the desired genes
  genes_features <- BiocGenerics:::subset(genes_assembly,
                                          symbol %in% genes)

  # Checking genes that could not be found and returned an error message
  `%!in%` <- base::Negate(`%in%`)

  all_genes <- genes_assembly$symbol %>%
    unlist() %>%
    unname()

  missing_genes <- genes[genes %!in% all_genes]

  if (!rlang::is_empty(missing_genes)) {
    warning(
      base::paste(
        "Genes:",
        paste(missing_genes,
              collapse = ", "),
        ",could not be found. Maybe you need to use a different gene alias?"
      )
    )
  }

  #finding overlaps
  olaps <-
    suppressWarnings(GenomicRanges::findOverlaps(genes_features,
                                                 ranges,
                                                 ignore.strand = TRUE))

  # creating a data frame that will contain the genes and positions (index) in the
  # pipeline ranges.
  # some genes might overlap more than one range (more than one bin), in this case
  # only one will be kept
  df <-
    tibble::tibble(gene = as.character(genes_features$symbol[S4Vectors::queryHits(olaps)]),
                   pos = S4Vectors::subjectHits(olaps)) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  # checking for genes that might have been blacklisted from the varbin pipeline
  blk_list <- genes[genes %!in% missing_genes]
  blk_list <- blk_list[blk_list %!in% df$gene]

  if (!rlang::is_empty(blk_list)) {
    warning(
      base::paste(
        "Genes:",
        paste(blk_list,
              collapse = ", "),
        "are in blacklisted regions of the Varbin pipeline and can't be plotted."
      )
    )
  }


  # obtaining seg ratios and sbsetting for the genes

  if (assay == 'segment_ratios') {
    seg_data <- segment_ratios(scCNA)
  }

  if (assay == 'integer') {
    seg_data <- SummarizedExperiment::assay(scCNA, 'integer')
  }

  seg_data_genes <- seg_data[df$pos, ] %>%
    dplyr::mutate(gene = df$gene)

  #long format for plotting
  seg_long <- tidyr::gather(seg_data_genes,
                            key = "sample",
                            value = "segratio", -gene)
  #plotting

  p <-
    ggplot2::ggplot(seg_long, aes(
      x = forcats::fct_reorder(gene, segratio),
      y = segratio
    )) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    my_theme

  # geom barplot plots a frequency barplot of each copy number state
  if (geom == 'barplot') {

     p <- ggplot2::ggplot(seg_long, aes(
        x = forcats::fct_reorder(gene, segratio),
        y = segratio
      )) +
      ggplot2::theme_classic() +
     geom_bar(aes(fill = as.factor(segratio)),
              position = 'fill', stat = 'identity') +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         breaks = scales::pretty_breaks(n = 10)) +
      scale_fill_viridis_d(option = 'inferno') +
      my_theme +
      ylab("percentage")

    print(p)
  }

  #geom violin
  if (geom == "violin" & is.null(label)) {
    p <- p +
      ggplot2::geom_violin()

    print(p)

  } else if (geom == "violin" & !is.null(label)) {
    warning("Coloring by label argument is only available for geom = 'swarm'.")

    p <- p +
      ggplot2::geom_violin()

    print(p)

  }

  # geom swarm
  if (geom == "swarm" & is.null(label)) {
    p <- p +
      ggbeeswarm::geom_quasirandom()

    print(p)

  } else if (geom == "swarm" & !is.null(label)) {
    # retrieving data

    lab <- dplyr::pull(metadata,
                       var = label)


    if (label == "superclones") {
      # coloring for discrete variable label
      p <- p +
        ggbeeswarm::geom_quasirandom(aes(fill = rep(lab,
                                                    each = length(df$gene))),
                                     shape = 21,
                                     size = 2.2,
                                     stroke = 0.2)

      color_lab <-
        list(ggplot2::scale_fill_manual(values = superclones_pal(),
                                        limits = force))

      p <- p + color_lab

      print(p)

    } else if (label == "subclones") {
      # coloring for discrete variable label
      p <- p +
        ggbeeswarm::geom_quasirandom(aes(fill = as.character(rep(
          lab,
          each = length(df$gene)
        ))),
        shape = 21,
        size = 2.2,
        stroke = 0.2)

      color_lab <-
        list(ggplot2::scale_fill_manual(values = subclones_pal(),
                                        limits = force))

      p <- p + color_lab

      print(p)

    } else if (is.numeric(lab))  {
      p <- p +
        ggbeeswarm::geom_quasirandom(aes(fill = rep(lab,
                                                    each = length(df$gene))),
                                     shape = 21,
                                     size = 2.2,
                                     stroke = 0.2)

      color_lab <- list(ggplot2::scale_color_viridis_c())

      p <- p + color_lab

      print(p)

    } else {
      # coloring for discrete variable label
      p <- p +
        ggbeeswarm::geom_quasirandom(aes(fill = rep(lab,
                                                    each = length(df$gene))),
                                     shape = 21,
                                     size = 2.2,
                                     stroke = 0.2)

      color_lab <- list(ggplot2::scale_fill_viridis_d(limits = force))

      p <- p + color_lab

      print(p)

    }

  }

}
