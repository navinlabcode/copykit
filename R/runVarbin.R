#' Read BAM file and perform varbin analysis
#'
#' runVarbin performs the variable binning (VarBin) algorithm to a series of BAM
#' files resulting from short-read sequencing.
#'
#' @author Darlan Conterno Minussi
#'
#' @param dir A path containing .BAM files from short-read sequencing.
#' @param genome A character indicating the choice of genome assembly.
#' @param resolution A character indicating the resolution for the scaffold
#'  of the VarBin method, i. e. the bin resulting bin size.
#' @param remove_Y A boolean when set to TRUE, removes information from the chrY
#' from the dataset.
#' @param vst A character indicating the variance stabilization transformation
#'  to be performed. See \link{runVst} details.
#' @param method A character indicating the segmentation method.
#' @param min_bincount A numerical indicating the minimum mean bin counts a
#' cell should have to remain in the dataset.
#' @param is_paired_end A boolean indicating if bam files are from single-read
#' or pair end sequencing.
#' @param seed A numeric scalar that sets the seed for CBS segmentation permutation
#' reproducibility.
#' @param alpha A numeric with the. significance levels for the test to accept
#' change-points for CBS segmentation. See \code{\link[DNAcopy]{segment}}.
#' @param gamma A numeric passed on to 'multipcf' segmentation. Penalty for each
#'  discontinuity in the curve, default is 40. See \code{\link[copynumber]{multipcf}}.
#' @param name A character with the name for the slot returned by \code{runVarbin}
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @return An scCNA object containing the bin counts, the ratios and the segment
#' ratios.
#'
#' @details runVarbin is a convenient wrapper for CopyKit's pre-processing module.
#' It runs \code{runCountReads}, \code{runVst} and, \code{runSegmentation}.
#'
#' \code{runCountReads} takes as input duplicate marked BAM files from whole
#' genome sequencing and runs the variable binning pipeline algorithm. Briefly,
#' the genome is split into pre-determined bins. The bin size is controlled by
#' the argument \code{bin_size}. By using VarBin, for a diploid cell, each bin
#' will receive equal amount of reads, controlling for mappability.
#' A lowess function is applied to perform GC correction across the bins. The argument
#' \code{genome} can be set to 'hg38' or 'hg19' to select the scaffolds genome
#' assembly.
#' Information regarding the alignment of the reads to the bins and from the bam
#' files are stored in the #' \code{\link[SummarizedExperiment]{colData}}.
#'
#' \code{runVst} performs variance stabilization to reduce the overdispersion
#' from the negative binomial distribution nature of the bin counts and reduce
#' technical bias. The argument \code{vst} controls the choice of the transformation
#' allowing either the Freeman-Tukey transformation by using the option 'ft' (recommended)
#' or a logarithmic transformation with the option 'log'. Using a 'log' transformation
#' may result in long segmentation times for a few cells with large breakpoint counts.
#'
#' \code{runSegmentation} Fits a piece-wise constant function to the transformed
#' the smoothed bin counts. Bin counts are smoothed with
#' \code{\link[DNAcopy]{smooth.CNA}}. Segmentation can be chosen to one of the
#' following:
#'
#' \itemize{
#'
#'    \item{CBS:} \code{runSegmentation} Fits a piece-wise constant function
#'    to the transformed the smoothed bin counts. Bin counts are smoothed with
#'    \code{\link[DNAcopy]{smooth.CNA}} using the Circular Binary Segmentation
#'    (CBS) algorithm from \code{\link[DNAcopy]{segment}} with default it applies
#'    undo.prune with value of 0.05.
#'
#'    \item{multipcf:} Performs the joint segmentation from the \code{copynumber}
#'    package using the \code{\link[copynumber]{multipcf}} function. By fitting
#'    piecewise constant curves with common breakpoints for all samples.
#'
#'
#' }
#'
#' The resulting segment means are further refined with MergeLevels to join
#' adjacent segments with non-significant differences in segmented means.
#'
#' @references
#' Navin, N., Kendall, J., Troge, J. et al. Tumour evolution inferred by single-cell
#' sequencing. Nature 472, 90–94 (2011). https://doi.org/10.1038/nature09807
#'
#' Baslan, T., Kendall, J., Ward, B., et al (2015). Optimizing sparse sequencing
#' of single cells for highly multiplex copy number profiling.
#' Genome research, 25(5), 714–724. https://doi.org/10.1101/gr.188060.114
#'
#' Olshen AB, Venkatraman ES, Lucito R, Wigler M. Circular binary segmentation
#' for the analysis of array-based DNA copy number data. Biostatistics.
#' 2004 Oct;5(4):557-72. doi: 10.1093/biostatistics/kxh008. PMID: 15475419.
#'
#' Freeman, M. F.; Tukey, J. W. (1950), "Transformations related to the angular
#' and the square root", The Annals of Mathematical Statistics,
#' 21 (4), pp. 607–611, doi:10.1214/aoms/1177729756, JSTOR 2236611
#'
#' Fryzlewicz, P. (2014). WILD BINARY SEGMENTATION FOR MULTIPLE CHANGE-POINT
#' DETECTION. The Annals of Statistics, 42(6), 2243-2281. Retrieved July 30,
#' 2021, from http://www.jstor.org/stable/43556493
#'
#' Nilsen G, Liestol K, Van Loo P, Vollan H, Eide M, Rueda O, Chin S, Russell R,
#' Baumbusch L, Caldas C, Borresen-Dale A, Lingjaerde O (2012). “Copynumber:
#' Efficient algorithms for single- and multi-track copy number segmentation.”
#' BMC Genomics, 13(1), 591.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' copykit_obj <- runVarbin("~/path/to/bam/files/", remove_Y = TRUE)
#' }
#'
runVarbin <- function(dir,
                      genome = c("hg38", "hg19"),
                      resolution = c(
                          "200kb",
                          "50kb",
                          "100kb",
                          "175kb",
                          "250kb",
                          "500kb",
                          "1Mb",
                          "2.5Mb"
                      ),
                      remove_Y = FALSE,
                      is_paired_end = FALSE,
                      method = c("CBS", "multipcf"),
                      vst = c("ft", "log"),
                      seed = 17,
                      min_bincount = 10,
                      alpha = 1e-9,
                      gamma = 40,
                      name = "segment_ratios",
                      BPPARAM = bpparam()) {
    genome <- match.arg(genome)
    vst <- match.arg(vst)
    method <- match.arg(method)
    resolution <- match.arg(resolution)

    copykit_object <- runCountReads(dir,
        genome = genome,
        resolution = resolution,
        remove_Y = remove_Y,
        min_bincount = min_bincount,
        is_paired_end = is_paired_end,
        BPPARAM = BPPARAM
    )

    copykit_object <- runVst(copykit_object,
        transformation = vst
    )

    copykit_object <- runSegmentation(copykit_object,
        seed = seed,
        name = name,
        alpha = alpha,
        gamma = gamma,
        method = method,
        BPPARAM = BPPARAM
    )

    copykit_object <- logNorm(copykit_object,
        transform = "log"
    )

    return(copykit_object)
}
