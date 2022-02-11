# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Documentation for internal data objects
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' hg38_grangeslist
#'
#' @name hg38_grangeslist
#' @aliases hg38_grangeslist
#' @docType data
#' @return Contains a GrangesList object with the scaffolds for each of the
#' resolutions used by runVarbin, runCountReads and runSegmentation on the
#' hg38 genome assembly.
#' @rdname data
#' @keywords internal
"hg38_grangeslist"

#' hg19_rg
#'
#' @name hg19_rg
#' @aliases hg19_rg
#' @docType data
#' @return Contains a GrangesList object with the scaffolds for each of the
#' resolutions used by runVarbin, runCountReads and runSegmentation on the hg19
#' genome assembly.
#' @rdname data
#' @keywords internal
"hg38_grangeslist"


#' copykit_obj_rle
#'
#' @name copykit_obj_rle
#' @aliases copykit_obj_rle
#' @docType data
#' @details Contains a Rle object with the segment ratios values of the copykit
#' BL1 example dataset. It is used by the functions copykit_example
#' @rdname data
#' @keywords internal
"copykit_obj_filt_rle"

#' copykit_obj_filt_rle
#'
#' @name copykit_obj_filt_rle
#' @aliases copykit_obj_filt_rle
#' @docType data
#' @details Contains a Rle object with the segment ratios values of the copykit
#' BL1 example filtered dataset. It is used by the functions
#' copykit_example_filtered()
#' @rdname data
#' @keywords internal
"copykit_obj_rle"

#' copykit_obj_filt_umap
#'
#' @name copykit_obj_filt_umap
#' @aliases copykit_obj_filt_umap
#' @docType data
#' @details Contains the umap reduced dimension for the BL1 dataset as generated
#' for the CopyKit manuscript
#' @rdname data
#' @keywords internal
"copykit_obj_filt_umap"

#' hg19_genes
#'
#' @name hg19_genes
#' @aliases hg19_genes
#' @docType data
#' @details Contains the GrangesObject for the genomic positions of genes in the
#' hg19 genome assembly
#' @source library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' @rdname data
#' @keywords internal
"hg19_genes"

#' hg38_genes
#'
#' @name hg38_genes
#' @aliases hg38_genes
#' @docType data
#' @details Contains the GrangesObject for the genomic positions of genes in the
#' hg38 genome assembly
#' @source library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' @rdname data
#' @keywords internal
"hg38_genes"
