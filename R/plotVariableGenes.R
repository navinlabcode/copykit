#' plotVariableGenes
#'
#' Visualization for the most variable genes found with \code{findVariableGenes}.
#'
#' @param scCNA scCNA object.
#' @param n A numeric defining how many variable genes will be plotted.
#'
#' @details \code{plotVariableGenes} plots the list of genes that was found
#' using \code{findVariableGenes}.
#'
#' @seealso \code{\link{findVariableGenes}}
#'
#' @importFrom BiocGenerics subset
#' @importFrom S4Vectors metadata subjectHits queryHits
#' @import ggplot2
#'
#' @return
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example()
#' copykit_obj <- findNormalCells(copykit_obj)
#' copykit_obj <- copykit_obj[,colData(copykit_obj)$is_normal == "FALSE"]
#' copykit_obj <- filterCells(copykit_obj)
#' copykit_obj <- copykit_obj[,colData(copykit_obj)$filtered == "kept"]
#' copykit_obj <- findVariableGenes(copykit_obj)
#' plotVariableGenes(copykit_obj)
#'
plotVariableGenes <- function(scCNA,
                              n = 30) {

  # checks
  if (is.null(S4Vectors::metadata(scCNA)$hvg)) {
    stop("Run findVariableGenes() first.")
  }

  # extracting data
  hvg_obj <- S4Vectors::metadata(scCNA)$hvg

  if (n > length(hvg_obj)) {
    stop('length of n must be smaller or equal to length of metadata(scCNA)$hvg.')
  }

  pca_df <- attr(hvg_obj, 'pca_df')
  hvg_obj <- hvg_obj[1:n]
  pca_df <- pca_df[hvg_obj,]

  pca_df <- pca_df %>%
    dplyr::mutate(gene = as.factor(gene)) %>%
    dplyr::mutate(gene = forcats::fct_reorder(gene, abs(pca_df$p1)))


  # theme setup
  my_theme <- list(
    ggplot2::theme(
      axis.title.x = element_text(colour = "gray28", size = 20),
      axis.text.x = element_text(size = 10),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(colour = "gray28", size = 20),
      axis.text.y = element_text(size = 10),
      # axis.line.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    )
  )

  p <- ggplot(pca_df, aes(x = gene, y = abs(p1))) +
    geom_segment(aes(
      x = gene,
      xend = gene,
      y = 0,
      yend = abs(p1)
    )) +
    geom_point(size = 4,
               fill = "#21908C",
               shape = 21) +
    theme_classic() +
    coord_flip()  +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, max(abs(pca_df$p1) + 0.02)),
      breaks = c(0, max(abs(pca_df$p1) - 0.05)),
      labels = c("Less variable", "More variable")
    ) +
    labs(x = "",
         y = "") +
    my_theme

  print(p)

}
