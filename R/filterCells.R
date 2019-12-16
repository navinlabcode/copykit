#' Filter noise cells
#'
#' filterCells uses a k-nearest-neighbor approach to remove cells with random CNA profiles, largely due to noise data. It calculates a correlation matrix and sets a resolution below which non neighbors will be classified as noise cells.
#' Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#'
#' @author Hua-Jun Wu
#'
#' @param scCNA scCNA object.
#' @param k K-nearest-neighbor, defaults to 5.
#' @param resolution Set's how strict the correlation cut off will be. Defaults to 0.8.
#'
#' @return Adds a filtered cells label to the scCNA metadata. Cells that pass the filtering criteria receive the label "kept", whereas cells that do not pass the filtering criteria receive the label "removed".
#' @export
#'
#' @examples
#'
#'

filterCells <- function(scCNA,
                        k = 5,
                        resolution = 0.8) {
  k = k
  seg <- copykit::segment_ratios(scCNA)

  message("Calculating correlation matrix.")
  dst = cor(seg)
  dst_knn_df = apply(as.matrix(dst), 1, function(x){
    mean(sort(x, decreasing=T)[2:(k+1)])
  }) %>%
    tibble::enframe(name = "sample",
                    value = "cor")

  dst_knn_df <- dst_knn_df %>%
    dplyr::mutate(filtered = dplyr::case_when(
      cor >= resolution ~ "kept",
      cor < resolution ~ "removed"
    ))

  message("Adding information to metadata. Access with SummarizedExperiment::colData(scCNA).")
  if (identical(colData(scCNA)$sample, names(dst_knn))) {
    colData(scCNA)$filtered <- dst_knn_df$filtered
  } else stop("Sample names do not match metadata sample info. Check colData(scCNA).")

  message("Done.")
  return(scCNA)

}

