#' runPca()
#'
#' Creates a pca embedding using the package uwot from the segment ratios values
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA The CopyKit object.
#' @param assay String with the name of the assay to pull data from to make the
#'embedding.
#' @param name String specifying the name to be used to store the result in the
#' reducedDims of the output.
#' @param scale A logical value indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place.
#' @param ... Additional parameters passed to \code{\link[stats]{prcomp}}.
#'
#' @importFrom stats prcomp
#' @importFrom SummarizedExperiment assay
#'
#' @return A reduced dimension representation with pca in the slot
#' \code{reducedDim} from scCNA object. Access reduced dimensions slot with:
#' \code{reducedDim(scCNA, 'PCA', withDimnames = FALSE)}
#' @export
#'
#' @examples
#' set.seed(1000)
#' copykit_obj <- copykit_example_filtered()[,sample(100)]
#' copykit_obj <- runPca(copykit_obj)

runPca <- function(scCNA,
                    assay = "logr",
                    name = "PCA",
                   scale = FALSE,
                    ...) {


  seg_data <- t(SummarizedExperiment::assay(scCNA, assay)) %>%
    as.data.frame()

  message(paste("Using assay:", assay))
  message(paste("Embedding data with PCA."))

  pca <- prcomp(seg_data,
                     scale. = scale,
                    ...)

  # Saving results
  dat_pca <- pca$x
  rownames(dat_pca) <- rownames(seg_data)
  var_explained <- pca$sdev^2
  attr(dat_pca, "var_explained") <- var_explained
  rownames(pca$rotation) <- colnames(seg_data)
  attr(dat_pca, "rotation") <- pca$rotation

  SingleCellExperiment::reducedDim(scCNA, type = name) <- dat_pca





  message(
    "Access reduced dimensions slot with: SingleCellExperiment::reducedDim(scCNA, 'pca')."
  )
  message("Done.")

  return(scCNA)

}
