#' Read BAM file and perform varbin analysis
#'
#' runVarbin performs the variable binning (VarBin) algorithm to a series of BAM files resulting from short-read sequencing.
#'
#' @author Darlan Conterno Minussi
#'
#' @param dir A path for the directory containing .BAM files from short-read sequencing.
#' @param genome Choice of genome assembly. Default: 'hg38'.
#' @param bin_size The resolution of the VarBin method. Default: '200kb'.
#' @param remove_Y (default == FALSE) If set to TRUE, removes information from the chrY from the dataset.
#' @param vst Character. Transformation to be performed. More details in \link{copykit::runVst}
#' @param method Character. Segmentation method of choice.
#' @param seed Numeric. Set seed for CBS segmentation permutation reproducibility.
#' @param slot Character. Target slot for the resulting segment ratios.
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details runVarbin is a convinient wrapper for CopyKit's pre-processing module.
#' It runs \link{runCountReads}, \link{runVst} and, \link{runSegmentation}.
#'
#' @return An scCNA object containing the bin counts, the ratios and the segment
#' ratios.
#'
#' @export
#'
#' @examples
#'
#'

runVarbin <- function(dir,
                      genome = "hg38",
                      bin_size = "200kb",
                      remove_Y = FALSE,
                      vst = 'ft',
                      seed = 17,
                      slot = 'segment_ratios',
                      BPPARAM = bpparam()) {


  copykit_object <- runCountReads(dir,
                                  bin_size = bin_size,
                                  remove_Y = remove_Y,
                                  BPPARAM = BPPARAM)

  copykit_object <- runVst(copykit_object,
                           transformation = vst)

  copykit_object <- runSegmentation(copykit_object,
                                    seed = seed,
                                    slot = slot,
                                    BPPARAM = BPPARAM)

  return(copykit_object)

}

