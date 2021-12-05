#' Run distance matrix calculations
#'
#' Performs distance matrix calculations that can be downstream used for
#' hierarchical clustering or phylogenetic analysis. Uses \code{amap::Dist()}
#' in order to parallelize distance calculations.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param metric distance metric passed to calculate the distance matrix (Defaults to "euclidean").
#' @param n_threads Number of threads used to calculate the distance matrix.
#' Passed to `amap::Dist`. As default it uses 1/4 of the detected cores available.
#'
#' @return A distance matrix in the slot \code{distMat} from scCNA object.
#' Access the distance matrix with: \code{copykit::distMat(scCNA, withDimnames = TRUE)}
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- runDistMat(copykit_obj)

runDistMat <- function(scCNA,
                       metric = "euclidean",
                       n_threads = 1) {
  # cores check
  if (n_threads < 1) {
    n_threads <- 1
  }

  message(paste("Calculating distance matrix with metric:", metric))
  message(paste("Using", n_threads, "cores."))

  seg_data <- t(segment_ratios(scCNA)) %>%
    as.data.frame()

  dist_mat <-
    amap::Dist(seg_data,
               method = metric,
               nbproc = n_threads)

  distMat(scCNA) <- dist_mat

  message("Access distance matrix with copykit::distMat()")
  message("Done.")

  return(scCNA)

}
