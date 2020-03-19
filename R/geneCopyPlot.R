#' Plot gene segment ratios
#'
#' geneCopyPlot is a way to visualize the number of copies for different genes across all the cells. It outputs a violin plot with the desired genes.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param genes Vector containing the HUGO Symbol for the genes of interest.
#' @param geom Which geom should be used for plotting, options are "violin" or "swarm". Defaults to "swarm"
#' @param label Color by an element of metadata. Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#'
#' @return A violin plot with the segment ratios for the genes of interest.
#' @importFrom BiocGenerics subset
#' @importFrom BiocGenerics `%in%`
#' @importFrom S4Vectors subjectHits
#' @importFrom forcats fct_reorder
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'

geneCopyPlot <- function(scCNA,
                         genes,
                         geom = "swarm",
                         label = NULL) {
  # checks
  # check if label exists
  if (!is.null(label)) {
    metadata <- as.data.frame(SummarizedExperiment::colData(scCNA))

    message(paste0("Coloring by: ", label))
  }

  if (!is.null(label) && !(label %in% colnames(metadata))) {
    stop(paste0("Label ", label, " is not a column of the scCNA object."))
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
    ylab("Segment ratio")
  )

  # getting ranges from scCNA object
  ranges <- SummarizedExperiment::rowRanges(scCNA)

  # subsetting to only the desired genes
  # hg19_genes genomic range is saved inside sysdata.rda to avoid loading a lot of annotation packages
  hg19_genes_features <- BiocGenerics:::subset(hg19_genes,
                                               SYMBOL %in% genes)

  # Checking genes that could not be found and returned an error message
  `%!in%` <- base::Negate(`%in%`)

  all_genes <- hg19_genes$SYMBOL %>%
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
    suppressWarnings(GenomicRanges::findOverlaps(hg19_genes_features,
                                                 ranges,
                                                 ignore.strand = TRUE))

  # creating a data frame that will contain the genes and positions (index) in the
  # pipeline ranges.
  # some genes might overlap more than one range (more than one bin), in this case
  # only one will be kept
  df <-
    tibble::tibble(
      gene = as.character(hg19_genes_features$SYMBOL[S4Vectors::queryHits(olaps)]),
      pos = S4Vectors::subjectHits(olaps),
    ) %>%
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
  seg_data <- segment_ratios(scCNA)

  seg_data_genes <- seg_data[df$pos,] %>%
    dplyr::mutate(gene = df$gene)

  #long format for plotting
  seg_long <- tidyr::gather(seg_data_genes,
                            key = "sample",
                            value = "segratio",-gene)
  #plotting

  p <-
    ggplot2::ggplot(seg_long, aes(
      x = forcats::fct_reorder(gene, segratio),
      y = segratio + 1e-3
    )) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    my_theme

    #geom violin

  if (geom == "violin" & is.null(label)) {
    p <- p +
      ggplot2::geom_violin()

  } else if (geom == "violin" & !is.null(label)) {
    warning("Coloring by label argument is only available for geom = 'swarm'.")

    p <- p +
      ggplot2::geom_violin()

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


   if (label == "major_clusters") {
      # coloring for discrete variable label
      p <- p +
        ggbeeswarm::geom_quasirandom(aes(color = rep(lab,
                                                     each = length(df$gene))))

      color_lab <-
        list(ggplot2::scale_color_manual(values = major_palette))

      p <- p + color_lab

      print(p)

    } else if (label == "minor_clusters") {
      # coloring for discrete variable label
      p <- p +
        ggbeeswarm::geom_quasirandom(aes(color = as.character(rep(lab,
                                                     each = length(df$gene)))))

      color_lab <-
        list(ggplot2::scale_color_manual(values = minor_palette))

      p <- p + color_lab

      print(p)

    } else if (is.numeric(lab))  {

        p <- p +
          ggbeeswarm::geom_quasirandom(aes(color = rep(lab,
                                                       each = length(df$gene))))

        color_lab <- list(ggplot2::scale_color_viridis_c())

        p <- p + color_lab

        print(p)

    } else {
      # coloring for discrete variable label
      p <- p +
        ggbeeswarm::geom_quasirandom(aes(color = rep(lab,
                                                     each = length(df$gene))))

      color_lab <- list(ggplot2::scale_color_viridis_d())

      p <- p + color_lab

      print(p)

    }

  }

}
