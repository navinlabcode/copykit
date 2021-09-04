#' Read BAM file and perform varbin analysis
#'
#' runVarbin performs the variable binning (VarBin) algorithm to a series of BAM
#' files resulting from short-read sequencing.
#'
#' @author Darlan Conterno Minussi
#'
#' @param dir A path for the directory containing .BAM files from short-read sequencing.
#' @param genome A character indicating the choice of genome assembly.
#' @param bin_size A character indicating the resolution desired for the scaffold
#'  of the VarBin method, i. e. the bin resulting bin size.
#' @param remove_Y A boolean when set to TRUE, removes information from the chrY
#' from the dataset.
#' @param vst A character indicating the variance stabilization transformation
#'  to be performed. See \link{runVst} details.
#' @param method A character indicating the segmentation method.
#' @param min_bincount A numerical indicating the minimum mean bin counts a
#' cell should have to remain in the dataset.
#' @param seed A numeric scalar that sets the seed for CBS segmentation permutation
#' reproducibility.
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
#' the smoothed bin counts. Bin counts are smoothed with \code{\link[DNAcopy]{smooth.CNA}}
#' using the Circular Binary Segmentation (CBS) algorithm from
#' \code{\link[DNAcopy]{segment}} with default it applies undo.prune with value of 0.05.
#' Or with Wild Binary Segmentation (WBS) from \code{\link[wbs]{wbs}}.
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
#' @export
#'
#' @examples
#'
#'

runVarbin <- function(dir,
                      genome = c("hg38", "hg19"),
                      bin_size = "200kb",
                      remove_Y = FALSE,
                      method = c('CBS', 'WBS', 'multipcf'),
                      vst = c("ft", "log"),
                      seed = 17,
                      min_bincount = 10,
                      name = 'segment_ratios',
                      BPPARAM = bpparam()) {

  genome <- match.arg(genome)
  vst <- match.arg(vst)
  method <- match.arg(method)

  copykit_object <- runCountReads(dir,
                                  genome = genome,
                                  bin_size = bin_size,
                                  remove_Y = remove_Y,
                                  min_bincount = min_bincount,
                                  BPPARAM = BPPARAM)

  copykit_object <- runVst(copykit_object,
                           transformation = vst)

  copykit_object <- runSegmentation(copykit_object,
                                    seed = seed,
                                    name = name,
                                    method = method,
                                    BPPARAM = BPPARAM)

  return(copykit_object)

}

