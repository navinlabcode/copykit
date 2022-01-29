#' plotScree
#'
#' Plots the variance explained by the different principal components
#'
#' @param scCNA The CopyKit object
#' @param ncomponents Number of principal components to plot.
#'
#' @return A ggplot object with The variance explained per principal component.
#' @export
#'
#' @importFrom scales percent_format
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- runPca(copykit_obj)
#' plotScree(copykit_obj)
#'
plotScree <- function(scCNA,
                      ncomponents = 20) {

  # sdev attribute is saved with the PCA redDim
  sdev <- attr(reducedDim(scCNA, "PCA"), 'var_explained')

  # Calculating variance explained
  ve <- sdev / sum(sdev)

  # Data frame for plotting
  df <- data.frame(pcacomponents = 1:ncomponents,
                   varexp = ve[1:ncomponents])

  p <- ggplot(df, aes(x = pcacomponents,
                      y = varexp)) +
    geom_point() +
    theme_classic() +
    theme(axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 16)) +
    labs(x = "PCA components",
         y = 'variance explained') +
    scale_y_continuous(labels = scales::percent_format())

  p

}
