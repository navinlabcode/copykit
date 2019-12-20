### all clases for the copykit package

###################################################################
# defining scCNA class
#'
#'  @export
#'  @import methods SingleCellExperiment
#'  @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SingleCellExperiment
#'  @importClassesFrom S4Vectors DataFrame SimpleList
#'

#' @export
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SingleCellExperiment
.scCNA <- setClass("scCNA", contains = "SingleCellExperiment")

#' @export
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SingleCellExperiment
scCNA <- function(segment_ratios,
                  ratios,
                  bin_counts,
                  ...) {
  cna <-
    SingleCellExperiment::SingleCellExperiment(list(segment_ratios = segment_ratios,
                                                    ratios = ratios,
                                                    bin_counts = bin_counts), ...)
  .scCNA(cna)
}
