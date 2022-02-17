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
#' @return A string vector with the HUGO genes in decreasing order of importance
#' stored to the \code{\link[S4Vectors]{metadata}}.
#'
#' @details \code{findVariableGenes} Runs \code{\link[stats]{prcomp}} to the
#' copy number states of the genes from the provided gene list and returns
#' the one that have the largest absolute variance as assesed by the
#' loadings of the first principal component.
#'
#' The resulting list of genes is stored within the metadata of the scCNA
#' object and can be accessed with \code{\link[S4Vectors]{metadata}}.
#'
#'
#' @importFrom SummarizedExperiment rowRanges
#' @importFrom BiocGenerics subset
#' @importFrom S4Vectors metadata subjectHits queryHits
#' @importFrom stats prcomp
#' @importFrom rlang is_empty
#' @importFrom GenomicRanges findOverlaps
#'
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- findVariableGenes(copykit_obj,
#'     genes = c("FHIT", "PTEN", "FOXO1", "BRCA1")
#' )
findVariableGenes <- function(scCNA,
                              genes,
                              assay = "logr",
                              top_n = 50) {
    # checks
    if (top_n > length(genes)) {
        top_n <- length(genes)
    }

    # obtaining df with genes positions
    # find_scaffold_genes in internals.R
    df <- find_scaffold_genes(scCNA,
        genes = genes
    )

    # obtaining data and subsetting
    seg_data <- SummarizedExperiment::assay(scCNA, assay)
    seg_data_genes <- seg_data[df$pos, ]

    rownames(seg_data_genes) <- df$gene

    # running principal component analysis
    pca_obj <- prcomp(t(seg_data_genes))

    pca_df <- data.frame(
        gene = names(pca_obj$rotation[, 1]),
        p1 = pca_obj$rotation[, 1]
    )

    pca_df <- pca_df[order(abs(pca_df$p1), decreasing = TRUE), ]

    hvg <- as.character(pca_df$gene[1:top_n])

    attr(hvg, "pca_pc1_loading") <- pca_obj$rotation[, 1]
    attr(hvg, "pca_df") <- pca_df

    S4Vectors::metadata(scCNA)$hvg <- hvg

    return(scCNA)
}
