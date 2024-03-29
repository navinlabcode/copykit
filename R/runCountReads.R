#' Aligns the reads from the BAM file to the variable binning pipeline.
#'
#' runCountReads performs the variable binning (VarBin) algorithm to a series of
#' BAM files resulting from short-read sequencing.
#'
#' @author Darlan Conterno Minussi
#'
#' @param dir A path for the directory containing BAM files from short-read
#' sequencing.
#' @param genome Name of the genome assembly. Default: 'hg38'.
#' @param resolution The resolution of the VarBin method. Default: '220kb'.
#' @param remove_Y (default == FALSE) If set to TRUE, removes information from
#' the chrY from the dataset.
#' @param min_bincount A numerical indicating the minimum mean bin counts a
#' cell should have to remain in the dataset.
#' @param is_paired_end A boolean indicating if bam files are from single-read
#' or pair end sequencing.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details \code{runCountReads} takes as input duplicate marked BAM files from
#' whole genome sequencing and runs the variable binning pipeline algorithm.
#' It is important that BAM files are duplicate marked. Briefly, the genome is
#' split into pre-determined bins. The bin size is controlled by the argument
#' \code{resolution}. By using VarBin, for a diploid cell, each bin will
#' receive equal amount of reads, controlling for mappability.
#' A lowess function is applied to perform GC correction across the bins.
#' The argument \code{genome} can be set to 'hg38' or 'hg19' to select the
#' scaffolds genome assembly. The scaffolds are GenomicRanges objects
#' Information regarding the alignment of the reads to the bins and from the bam
#' files are stored in the #' \code{\link[SummarizedExperiment]{colData}}.
#' \code{min_bincount} Indicates the minimum mean bincount a cell must present
#' to be kept in the dataset. Cells with low bincounts generally present bin
#' dropouts due to low read count that will be poorly segmented.
#'
#' @return A matrix of bin counts within the scCNA object that can be accessed
#' with \code{bincounts}
#'
#' #' @references
#' Navin, N., Kendall, J., Troge, J. et al. Tumour evolution inferred by
#' single-cell sequencing. Nature 472, 90–94 (2011).
#' https://doi.org/10.1038/nature09807
#'
#' Baslan, T., Kendall, J., Ward, B., et al (2015). Optimizing sparse sequencing
#' of single cells for highly multiplex copy number profiling.
#' Genome research, 25(5), 714–724. https://doi.org/10.1101/gr.188060.114
#'
#' @importFrom Rsubread featureCounts
#' @importFrom dplyr rename mutate relocate
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors DataFrame metadata
#'
#' @export
#'
#' @examples
#' \dontrun{
#' copykit_obj <- runCountReads("/PATH/TO/BAM/FILES")
#' }
#'
runCountReads <- function(dir,
                          genome = c("hg38", "hg19"),
                          resolution = c(
                              "220kb",
                              "55kb",
                              "110kb",
                              "195kb",
                              "280kb",
                              "500kb",
                              "1Mb",
                              "2.8Mb"
                          ),
                          remove_Y = FALSE,
                          min_bincount = 10,
                          is_paired_end = FALSE,
                          BPPARAM = bpparam()) {
    genome <- match.arg(genome)
    resolution <- match.arg(resolution)

    # bindings for NSE and data
    Chr <- chr <- strand <- GeneID <- NULL
    reads_assigned_bins <- reads_duplicates <- reads_total <- NULL

    files <-
        list.files(dir,
            pattern = "*.bam",
            full.names = TRUE,
            ignore.case = TRUE
        )

    if (!any(grepl(".bam", files, ignore.case = TRUE))) {
      stop("Directory does not contain .bam files.")
    }

    # managing .bai files
    if (any(grepl(".bai", files, ignore.case = TRUE))) {
      files <- files[!grepl(".bai", files)]
    }

    files_names <- basename(gsub(pattern = ".bam", "", files))

    # stop if files variable length is 0
    if (length(files) == 0) {
        stop("No .bam files detected.")
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # genomic ranges (varbin scaffolds)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Reading hg38 VarBin ranges
    if (genome == "hg38") {
        hg38_grangeslist <- hg38_grangeslist

        hg38_rg <- switch(resolution,
            "55kb" = hg38_grangeslist[["hg38_50kb"]],
            "110kb" = hg38_grangeslist[["hg38_100kb"]],
            "195kb" = hg38_grangeslist[["hg38_175kb"]],
            "220kb" = hg38_grangeslist[["hg38_200kb"]],
            "280kb" = hg38_grangeslist[["hg38_250kb"]],
            "500kb" = hg38_grangeslist[["hg38_500kb"]],
            "1Mb" = hg38_grangeslist[["hg38_1Mb"]],
            "2.8Mb" = hg38_grangeslist[["hg38_2Mb"]]
        )

        hg38_rg <- as.data.frame(hg38_rg)

        rg <- hg38_rg %>%
            dplyr::rename(chr = "seqnames") %>%
            dplyr::mutate(GeneID = 1:nrow(hg38_rg))

        if (remove_Y == TRUE) {
            rg <- dplyr::filter(
                rg,
                chr != "chrY"
            )
        }
    }

    # reading hg19 varbin ranges
    if (genome == "hg19") {
        rg <- hg19_rg %>%
            dplyr::mutate(GeneID = 1:nrow(hg19_rg))

        if (remove_Y == TRUE) {
            rg <- dplyr::filter(
                rg,
                chr != "chrY"
            )
        }
    }

    message(
        "Counting reads for genome ",
        genome,
        " and resolution: ",
        resolution
    )

    varbin_counts_list_all_fields <-
        suppressMessages(
            BiocParallel::bplapply(
                files,
                Rsubread::featureCounts,
                ignoreDup = TRUE,
                countMultiMappingReads = FALSE,
                annot.ext = rg,
                useMetaFeatures = FALSE,
                verbose = FALSE,
                isPairedEnd = is_paired_end,
                BPPARAM = BPPARAM
            )
        )

    names(varbin_counts_list_all_fields) <- files_names

    varbin_counts_list <- lapply(
        varbin_counts_list_all_fields,
        "[[",
        1
    )

    varbin_counts_list <- lapply(
        varbin_counts_list,
        as.vector
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # filtering for minimal mean bin count
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # obtaining the index of the ones that FAIL to meet the min_bincount arg
    min_bc <- which(vapply(varbin_counts_list, mean, numeric(1)) < min_bincount)
    # subsetting counts list and main counts list

    if (length(min_bc) > 0) {
        varbin_counts_list <- varbin_counts_list[-min_bc]
        varbin_counts_list_all_fields <- varbin_counts_list_all_fields[-min_bc]
        message(
            length(min_bc), " bam files had less than ", min_bincount,
            " mean bincounts and were removed."
        )
    }

    # LOWESS GC normalization

    message("Performing GC correction.")

    varbin_counts_list_gccor <-
        BiocParallel::bplapply(varbin_counts_list, function(x) {
            gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
            gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
            exp(log(x) - gc_cor_z$y) * median(x)
        },
        BPPARAM = BPPARAM
        )

    varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)

    # filtering low read counts where the sum of bins does not reach more than 0
    good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])

    varbin_counts_df <- varbin_counts_df[good_cells]

    rg <- rg %>%
        dplyr::select(-strand, -GeneID)

    rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,
        ignore.strand = TRUE,
        keep.extra.columns = TRUE
    )


    cna_obj <- CopyKit(
        assays = list(bincounts = varbin_counts_df),
        rowRanges = rg_gr
    )

    # Adding genome and resolution information to metadata
    S4Vectors::metadata(cna_obj)$genome <- genome
    S4Vectors::metadata(cna_obj)$resolution <- resolution

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:01 2021
    # ADDING READS METRICS TO METADATA
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sun Feb 14 20:55:24 2021

    varbin_reads_list <- lapply(
        varbin_counts_list_all_fields,
        "[[",
        4
    )

    # saving info and removing columns from list elements
    metadata_info_names <- varbin_reads_list[[1]][c(1, 2, 8, 9, 12, 14), 1]
    metadata_info_names <-
        c(
            "reads_assigned_bins",
            "reads_unmapped",
            "reads_duplicates",
            "reads_multimapped",
            "reads_unassigned",
            "reads_ambiguous"
        )

    varbin_reads_info <-
        lapply(seq_along(varbin_reads_list), function(x) {
            # RSubread seems to change underlines to dot on some cases
            # Have to make more complicated lapply to extract the name of the list
            # and guarantee that the cell is properly named
            name <- names(varbin_reads_list)[[x]]
            df <- varbin_reads_list[[x]][c(1, 2, 8, 9, 12, 14), -1, drop = FALSE]
            names(df) <- name
            df
        })

    names(varbin_reads_list) <- names(varbin_counts_list_all_fields)

    bam_metrics <- dplyr::bind_cols(varbin_reads_info)

    # making sure metrics match varbin_counts_df
    bam_metrics <- bam_metrics[good_cells]
    rownames(bam_metrics) <- metadata_info_names
    bam_metrics <- as.data.frame(t(bam_metrics))

    # adding total
    reads_tot <- rowSums(bam_metrics)

    bam_metrics$sample <- rownames(bam_metrics)
    bam_metrics <-
        dplyr::relocate(bam_metrics, sample, .before = reads_assigned_bins)

    bam_metrics <- bam_metrics %>%
        dplyr::mutate(
            reads_total = reads_tot,
            percentage_duplicates = round(reads_duplicates / reads_total, 3)
        )

    if (sum(bam_metrics$reads_duplicates) == 0) {
        warning(
            "runCountReads did not detect any duplicate reads.
      Make sure your input bam files have duplicates marked.",
            call. = FALSE,
            noBreaks. = TRUE
        )
    }

    # adding to metadata
    SummarizedExperiment::colData(cna_obj) <-
        S4Vectors::DataFrame(bam_metrics)
    colnames(cna_obj) <- names(varbin_counts_df)

    return(cna_obj)
}
