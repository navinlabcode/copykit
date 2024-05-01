#' plotGeneCopy
#'
#' Visualization for gene-wise copy number states
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param genes A vector of strings containing the HUGO Symbol for the gene
#' of interest.
#' @param geom A string with the geom for plotting.
#' @param label A string with the name of the column from
#' \code{\link[SummarizedExperiment]{colData}} to color the points
#' @param facet A string with the name of the column from
#' \code{\link[SummarizedExperiment]{colData}} to separate the plot into facets.
#' @param dodge.width A numeric that adds dodge between the label categories.
#' @param assay String with the name of the assay for plotting.
#'
#' @details plotGeneCopy finds overlaps of the varbin scaffolds genomic ranges
#' which can be accessed with \code{\link[SummarizedExperiment]{rowRanges}}
#' with the genes from the human genome assemblies of either hg19 or hg38, 
#' or the mouse genome assembly mm10. All are from TxDb.Hsapiens.UCSC packages.
#' hg19 package \code{TxDb.Hsapiens.UCSC.hg19.knownGene}.
#' hg38 package \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#' mm10 package \code{TxDb.Mmusculus.UCSC.mm10.knownGene}.
#'
#' If the argument geom is set to 'barplot' plotGeneCopy calculates gene-wise
#'  frequencies of each copy number state for the selected genes across all of
#'  the cells. Geom 'barplot' can only be used with the argument
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
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- findClusters(copykit_obj)
#' plotGeneCopy(copykit_obj, genes = c("FHIT", "PTEN", "FOXO1", "BRCA1"))
#' plotGeneCopy(copykit_obj,
#'     genes = c("FHIT", "PTEN", "FOXO1", "BRCA1"),
#'     label = "subclones"
#' )
#' plotGeneCopy(copykit_obj,
#'     genes = c("FHIT", "PTEN", "FOXO1", "BRCA1"),
#'     label = "subclones", dodge.width = 0.8
#' )
plotGeneCopy <- function(scCNA,
                         genes,
                         geom = c("swarm", "barplot", "violin"),
                         label = NULL,
                         facet = NULL,
                         dodge.width = 0,
                         assay = "segment_ratios") {
    geom <- match.arg(geom)

    # bindings for NSE
    gene <- segratio <- plot_label <- plot_facet <- NULL

    # checks
    if (geom == "barplot" && assay != "integer") {
        stop("Argument geom 'barplot' can only be used with assay == 'integer'")
    }

    # check for duplicated column names
    if (any(duplicated(colnames(scCNA)))) {
        colnames(scCNA) <- make.names(colnames(scCNA),
                                            unique = TRUE)

        warning("Warning: Detected and corrected duplicated cell names.")

    }

    # check if label is provided
    if (!is.null(label)) {
        metadata <- as.data.frame(SummarizedExperiment::colData(scCNA))
        message("Coloring by: ", label)
        metadata$sample <- colnames(scCNA)
    }

    # check if label exists
    if (!is.null(label) && !(label %in% colnames(metadata))) {
        stop("Label ", label, " is not a column of the scCNA object.")
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
            legend.position = "right",
            legend.title = element_blank(),
            legend.text = element_text(size = 16)
        ),
        xlab(""),
        ylab("segment ratio")
    )

    # obtaining df with genes positions
    # find_scaffold_genes in internals.R
    df <- find_scaffold_genes(scCNA,
                              genes = genes
    )

    # obtaining seg ratios and sbsetting for the genes

    if (assay == "segment_ratios") {
        seg_data <- segment_ratios(scCNA)
    }

    if (assay == "integer") {
        seg_data <- SummarizedExperiment::assay(scCNA, "integer")
    }

    seg_data_genes <- seg_data[df$pos,]
    seg_data_genes$gene <- df$gene



    # long format for plotting
    seg_long <- tidyr::gather(seg_data_genes,
                              key = "sample",
                              value = "segratio", -gene
    )

    # adding metadata if provided
    if (!is.null(label)) {
        metadata_info <- metadata[c("sample", label)]
        names(metadata_info) <- c("sample", "plot_label")
        seg_long <- merge(seg_long, metadata_info)
    }

    # adding facet if provided
    if (!is.null(facet)) {
        metadata_info_facet <- metadata[c("sample", facet)]
        names(metadata_info_facet) <- c("sample", "plot_facet")
        seg_long <- merge(seg_long, metadata_info_facet)
    }

    # plotting

    p <-
        ggplot2::ggplot(seg_long, aes(
            x = forcats::fct_reorder(gene, segratio),
            y = segratio
        )) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
        my_theme

    # geom barplot plots a frequency barplot of each copy number state
    if (geom == "barplot") {
        p <- ggplot2::ggplot(seg_long, aes(
            x = forcats::fct_reorder(gene, segratio),
            y = segratio
        )) +
            ggplot2::theme_classic() +
            geom_bar(aes(fill = as.factor(segratio)),
                     position = "fill", stat = "identity"
            ) +
            scale_y_continuous(
                labels = scales::percent_format(accuracy = 1),
                breaks = scales::pretty_breaks(n = 10)
            ) +
            scale_fill_viridis_d(option = "inferno") +
            my_theme +
            ylab("percentage")

        # return plot
        return(p)
    }

    # geom violin
    if (geom == "violin" & is.null(label)) {
        p <- p +
            ggplot2::geom_violin()

        return(p)
    } else if (geom == "violin" & !is.null(label)) {
        warning("Coloring by label argument is only available for geom = 'swarm'.")

        p <- p +
            ggplot2::geom_violin()

        return(p)
    }

    # geom swarm
    if (geom == "swarm" & is.null(label)) {
        p <- p +
            ggbeeswarm::geom_quasirandom()

        p
    } else if (geom == "swarm" & !is.null(label)) {
        # retrieving data

        if (label == "superclones") {
            # coloring for discrete variable label
            p <- p +
                ggbeeswarm::geom_quasirandom(aes(fill = plot_label),
                                             shape = 21,
                                             size = 2.2,
                                             dodge.width = dodge.width,
                                             stroke = 0.2
                )

            color_lab <-
                list(ggplot2::scale_fill_manual(
                    values = superclones_pal(),
                    limits = force
                ))

            p <- p + color_lab
        } else if (label == "subclones") {
            # coloring for discrete variable label
            p <- p +
                ggbeeswarm::geom_quasirandom(aes(fill = plot_label),
                                             shape = 21,
                                             size = 2.2,
                                             dodge.width = dodge.width,
                                             stroke = 0.2
                )

            color_lab <-
                list(ggplot2::scale_fill_manual(
                    values = subclones_pal(),
                    limits = force
                ))

            p <- p + color_lab
        } else if (is.numeric(metadata[[label]])) {
            p <- p +
                ggbeeswarm::geom_quasirandom(aes(fill = plot_label),
                                             shape = 21,
                                             size = 2.2,
                                             stroke = 0.2
                )

            color_lab <- list(ggplot2::scale_color_viridis_c())

            p <- p + color_lab
        } else {
            # coloring for discrete variable label
            p <- p +
                ggbeeswarm::geom_quasirandom(aes(fill = plot_label),
                                             shape = 21,
                                             size = 2.2,
                                             dodge.width = dodge.width,
                                             stroke = 0.2
                )

            color_lab <- list(ggplot2::scale_fill_viridis_d(limits = force))

            p <- p + color_lab
        }

        # adding facets if provided
        if (!is.null(facet)) {
            p <- p + facet_wrap(vars(plot_facet))
        }

        p
    }
}
