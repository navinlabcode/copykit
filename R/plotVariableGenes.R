#' plotVariableGenes
#'
#' Uses PCA to detect the most variable genes.
#'
#' @param scCNA scCNA object.
#' @param genes A vector of strings containing the HUGO Symbol for the gene
#' of interest.
#' @param assay String with the name of the assay to pull data from to calculate
#' integers.
#' @param top_n A numeric that if provided returns the top N variable genes.
#'
#' @return
#' @export
#'
#' @examples
plotVariableGenes <- function(scCNA,
                              genes,
                              assay = "segment_ratios",
                              top_n = NULL) {

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
      axis.text.x = element_text(size = 15),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour = "gray28", size = 20),
      axis.text.y = element_text(size = 6),
      axis.line.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    )
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

  seg_data_genes <- seg_data[df$pos,]

  rownames(seg_data_genes) <- df$gene

  pca_obj <- prcomp(t(seg_data_genes))

  pca_df <- data.frame(gene = names(pca_obj$rotation[, 1]),
                       p1 = pca_obj$rotation[, 1])

  pca_df <- pca_df %>%
    dplyr::mutate(gene = as.factor(gene)) %>%
    dplyr::mutate(gene = forcats::fct_reorder(gene, abs(p1)))

  pca_df <- pca_df[order(abs(pca_df$p1), decreasing = TRUE),]

  p <- ggplot(pca_df, aes(x = gene, y = abs(p1))) +
    geom_segment(aes(
      x = gene,
      xend = gene,
      y = 0,
      yend = abs(p1)
    )) +
    geom_point(size = 5,
               fill = "#21908C",
               shape = 21) +
    theme_classic() +
    coord_flip()  +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, max(abs(pca_df$p1) + 0.02)),
      breaks = c(0, max(abs(pca_df$p1) - 0.1)),
      labels = c("Less variable", "More variable")
    ) +
    labs(x = "",
         y = "") +
    my_theme

  if (!is.null(top_n)) {
    return(pca_df$gene[1:top_n])
  }

  print(p)

}
