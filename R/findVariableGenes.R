#' findVariableGenes
#'
#' Find the most variable genes in the dataset.
#'
#' @param scCNA scCNA object.
#' @param genes A vector of strings containing the HUGO Symbol for the gene
#' of interest.
#' @param assay String with the name of the assay to pull data with the copy
#' number states for each gene.
#' @param top_n A numeric defining how many variable genes will be returned.
#'
#' @details \code{findVariableGenes} Runs \code{\link[stats]{prcomp}} to the
#' copy number states of the genes from the provided gene list and returns
#' the one that have the largest absolute variance as assesed by the
#' loadings of the first principal component.
#'
#' The resulting list of genes is stored within the metadata of the scCNA
#' object and can be accessed with \code{\link[S4Vectors]{metadata}}.
#'
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom BiocGenerics subset
#' @importFrom S4Vectors metadata subjectHits queryHits
#' @importFrom tibble tibble
#' @importFrom rlang is_empty
#' @importFrom GenomicRanges findOverlaps
#'
#' @return
#' @export
#'
#' @examples
findVariableGenes <- function(scCNA,
                              genes,
                              assay = "segment_ratios",
                              top_n = 50) {
  # checks
  if (top_n > length(genes)) {
    stop("top_n length must be smaller than genes argument length.")
  }

  # genome assembly
  if (S4Vectors::metadata(scCNA)$genome == "hg19") {
    genes_assembly <- hg19_genes
  }

  if (S4Vectors::metadata(scCNA)$genome == "hg38") {
    genes_assembly <- hg38_genes
  }

  # getting ranges from scCNA object
  ranges <- SummarizedExperiment::rowRanges(scCNA)

  # subsetting to only the desired genes
  genes_features <- BiocGenerics:::subset(genes_assembly,
                                          symbol %in% genes)

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

  # obtaining data and subsetting
  seg_data <- SummarizedExperiment::assay(scCNA, assay)
  seg_data_genes <- seg_data[df$pos, ]

  rownames(seg_data_genes) <- df$gene

  # running principal component analysis
  pca_obj <- prcomp(t(seg_data_genes))

  pca_df <- data.frame(gene = names(pca_obj$rotation[, 1]),
                       p1 = pca_obj$rotation[, 1])

  pca_df <- pca_df[order(abs(pca_df$p1), decreasing = TRUE), ]

  hvg <- as.character(pca_df$gene[1:top_n])

  attr(hvg, 'pca_pc1_loading') <- pca_obj$rotation[, 1]
  attr(hvg, 'pca_df') <- pca_df

  S4Vectors::metadata(scCNA)$hvg <- hvg

  return(scCNA)

}
