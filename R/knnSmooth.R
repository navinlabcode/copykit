#' knnSmooth
#'
#' Smooth bincounts based on k nearest neighbors.
#'
#' @author Darlan Conterno Minussi
#' @author Runmin Wei
#'
#' @param scCNA The CopyKit object.
#' @param k A numeric with the k nearest neighbor value for smoothing
#' @param BPPARAM  A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @return The CopyKit object with an assay smoothed_bincounts
#'
#' @details This function uses a k-nearest neighbors approach to smooth cells
#' raw bincounts. To do so, the k-nearest neighbors are calculated with
#' \code{\link[BiocNeighbors]{findKNN}}. The bincounts of the k-nearest neighbors
#' for each cell are tallied and an assay called smoothed_bincounts is added to
#' \code{\link{assay}}. After, \code{\link{runVst}} and
#' \code{\link{runSegmentation}}. Are re-run by \code{knnSmooth}.
#'
#' This function results in a trade-off for the elimination of noise at the cost
#' of risk of loss of subclonal structure. To minimize the risk of subclonal
#' structure loss we recommend using the very small values of k.
#'
#' This function should be followed by applying \code{\link{runVst}} and
#' \code{\link{runSegmentation}} to the CopyKit object.
#'
#' @importFrom BiocNeighbors findKNN
#'
#' @export
#'
#' @examples
#' copykit_obj <- mock_bincounts(ncells = 10)
#' copykit_obj <- runSegmentation(copykit_obj)
#' copykit_obj <- knnSmooth(copykit_obj)
#'
#'
knnSmooth <- function(scCNA,
                      k = 3,
                      BPPARAM = bpparam()) {
  # setup data
  bin <- bincounts(scCNA)
  seg <- segment_ratios(scCNA)

  # finding neighbors
  message("Finding neighbors.")
  neighbors <- BiocNeighbors::findKNN(t(seg), k = k)

  message(paste("Smoothing cells using k =", k))
  # collect neighbors and sum counts
  smoothed_bins_list <- bplapply(seq_along(bin), function(i) {
    cells_neighbors_df <- bin[c(i, neighbors$index[i,])]
    smoothed_cell <- rowSums(cells_neighbors_df)
    smoothed_cell
  })

  # re-adding names
  names(smoothed_bins_list) <- colnames(scCNA)

  smoothed_bins_df <- as.data.frame(do.call(cbind,
                              smoothed_bins_list))

  # adding knn smoothed bins to assay and re-running Vst
  assay(scCNA, 'smoothed_bincounts') <- smoothed_bins_df

  # re-running vst and segmentation
  scCNA <- runVst(scCNA, assay = 'smoothed_bincounts')
  scCNA <- runSegmentation(scCNA)

  return(scCNA)

}
