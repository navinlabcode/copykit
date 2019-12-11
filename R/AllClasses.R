### all clases for the copykit package

###################################################################
# defining scCNA class
#'
#'  @export
#'  @import methods
#'  @importClassesFrom SummarizedExperiment RangedSummarizedExperiment SingleCellExperiment
#'  @importClassesFrom S4Vectors DataFrame SimpleList
#'
#'

.scCNA <- setClass("scCNA", contains = "SingleCellExperiment")

scCNA <- function(segment_ratios, ...) {
  cna <-
    SingleCellExperiment::SingleCellExperiment(list(segment_ratios = segment_ratios), ...)
  .scCNA(cna)
}
