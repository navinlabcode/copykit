#' Plot UMAP embedding
#'
#' Plots UMAP embedding stored in \code{reducedDims()} slot.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object. Access reduced dimensions slot with: \code{SingleCellExperiment::reducedDim(scCNA, 'umap')}. \code{plotUmap} searches for cluster information in the \code{SingleCellExperiment::colData()} metadata and colors the clusters according to that information.
#'
#' @export
#'
#' @examples
#'

plotUmap <- function(scCNA) {
  # theme setup
  my_theme <- list(theme(
    axis.title.x = element_text(colour = "gray28", size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(colour = "gray28", size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  ),xlab("umap1"),ylab("umap2"))

  # obtaining data from reducedDim slot
  if (!is.null(SingleCellExperiment::reducedDim(breast_tumor))) {
    umap_df <- SingleCellExperiment::reducedDim(scCNA, 'umap') %>%
      as.data.frame()

  } else
    stop("Reduced dimensions slot is null. Use runUmap() to create it.")

  if (is.null(SummarizedExperiment::colData(scCNA, "cluster"))) {
    message("No cluster information detected, use findClusters() to create it.")
    message("Plotting Umap.")

    ggplot(umap_df) +
      geom_point(aes(V1, V2)) +
      theme_classic() +
      my_theme

  } else {
    message("Plotting Umap.")
    message("Using colData(scCNA, 'cluster') cluster information.")

    ggplot(umap_df) +
      geom_point(aes(
        x = V1,
        y = V2,
        color = SummarizedExperiment::colData(scCNA, "cluster")
      )) +
      scale_color_manual(values = color_palette)
    theme_classic() +
      my_theme

  }


}