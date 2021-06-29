#' Plot UMAP embedding
#'
#' Plots UMAP embedding stored in \code{reducedDims()} slot.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param label Color by an element of metadata. Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#'
#' @return A reduced dimension representation with UMAP in the slot \code{reducedDim} from scCNA object. Access reduced dimensions slot with: \code{SummarizedExperiment::reducedDim(scCNA, 'umap', withDimnames = FALSE)}. \code{plotUmap} searches for cluster information in the \code{SummarizedExperiment::colData()} metadata and colors the clusters according to that information.
#'
#' @export
#'
#' @importFrom ggnewscale new_scale_color
#' @examples
#'

plotUmap <- function(scCNA,
                     label = NULL) {

  # retrieving metadata
  df <- as.data.frame(SummarizedExperiment::colData(scCNA))

  # theme setup
  my_theme <- list(
    ggplot2::theme(
      axis.title.x = element_text(colour = "gray28", size = 20),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(colour = "gray28", size = 20),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 14)
    ),
    xlab("umap1"),
    ylab("umap2")
  )

  # check if label exists
  if (!is.null(label)) {
    message(paste0("Coloring by: ", label))
  }

  if (!is.null(label) && !(label %in% colnames(df))) {
    stop(paste0("Label ", label, " is not a column of the scCNA object."))
  }

  # obtaining data from reducedDim slot
  if (!is.null(SingleCellExperiment::reducedDim(scCNA))) {
    umap_df <- SingleCellExperiment::reducedDim(scCNA, 'umap') %>%
      as.data.frame()

  } else
    stop("Reduced dimensions slot is null. Use runUmap() to create it.")

  message("Plotting Umap.")

  if (is.null(label) && is.null(SummarizedExperiment::colData(scCNA)$subclones)) {
    # if label is not provided and clusters were not run

    message("No cluster information detected, use findClusters() to create it.")

    ggplot(umap_df) +
      geom_point(aes(V1, V2)) +
      theme_classic() +
      my_theme

  } else if (is.null(label) && is.null(SummarizedExperiment::colData(scCNA)$superclones)) {
    # if label is not provided but findClusters was run for subclones

    message("Using colData(scCNA) subclones information.")

    ggplot(umap_df) +
      geom_point(aes(
        x = V1,
        y = V2,
        fill = as.factor(SummarizedExperiment::colData(scCNA)$subclones)
      ),
      size = 1.8,
      shape = 21) +
    scale_fill_manual(values = subclones_pal(),
                       name = "subclones",
                      limits = force) +
    theme_classic() +
      my_theme

  }  else if (is.null(label) && !is.null(SummarizedExperiment::colData(scCNA)$superclones)) {
    # if label is not provided but findClusters was run for superclones and subclones

    message("Using colData(scCNA) cluster information.")

    ggplot(umap_df) +
      geom_point(
        aes(
          x = V1,
          y = V2,
          color = SummarizedExperiment::colData(scCNA)$superclones
        ),
        alpha = 1,
        size = 5
      ) +
      scale_color_manual(values = superclones_pal(),
                         name = "superclones") +
      ggnewscale::new_scale_color() +
      geom_point(aes(
        x = V1,
        y = V2,
        fill = as.factor(SummarizedExperiment::colData(scCNA)$subclones)
      ),
      size = 1.8,
      shape = 21) +
      scale_fill_manual(values = subclones_pal(),
                        name = "subclones",
                        limits = force) +
      theme_classic() +
      my_theme

  } else if (!is.null(label)) {
    # if label is provided, coloring by the label

    lab <- dplyr::pull(df,
                       var = label)

        p <- ggplot(umap_df) +
      geom_point(aes(V1, V2,
                     fill = lab),
                 size = 1.8,
                 shape = 21) +
      theme_classic() +
      labs(fill = label) +
      my_theme

    # coloring by continuos variable
    if (is.numeric(lab)) {

      color_lab <- list(ggplot2::scale_fill_viridis_c())

      p <- p + color_lab

      print(p)

    } else {

      print(p)
    }

  }


}
