#' plotSuggestedK
#'
#' Uses the information from \code{\link{findSuggestedK}} to plot the values
#' of jaccard similarity from the tested k range on \code{\link{findSuggestedK}}.
#'
#' @param scCNA The scCNA object.
#' @param geom A character with the geom to be used for plotting.
#'
#' @details \code{\link{plotSuggestedK}} access the \code{\link[S4Vectors]{metadata}}
#' element suggestedK_df that is saved to the scDNA object after running
#' \code{\link{findSuggestedK}}. The dataframe is used for plotting either a
#' heatmap, when the argument geom = 'tile', or a dotplot when argument geom =
#' 'dotplot' or a boxplot when geom = 'boxplot'.
#'
#' \itemize{
#' \item{geom = 'boxplot':} Plots a boxplot of the jaccard similarities across
#' all clusters detected in the grid search. The large red points represent
#' the mean jaccard similarity.
#'
#' #' \item{geom = 'tile':} Plots a heatmap of the jaccard similarities across
#' all clusters detected in the grid search. The filling collors represent the
#' jaccard similarity value. Rows represent clusters and columns the k value.
#'
#' #' \item{geom = 'dotplot':} Plots a dotplot of the jaccard similarities across
#' all clusters detected in the grid search. Where the size of the dots represent
#' the jaccard similarity for each assesed k value.
#'
#' #' \item{geom = 'scatterplot':} Plots a scatterplot of the jaccard similarity
#' explained by the number of cells. Points are colored by subclone and lines
#' represent a linear regression across the points.
#'
#' }
#'
#' @return A ggplot2 object with the plot of different tested k values and their
#' jaccard similarity for each subclone
#'
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr mutate
#' @importFrom tidyr complete
#' @importFrom S4Vectors metadata
#' @importFrom gtools mixedsort
#'
#' @examples
plotSuggestedK <- function(scCNA,
                           geom = c('boxplot', 'tile', 'dotplot', 'scatterplot')) {

  geom <- match.arg(geom)

  df <- S4Vectors::metadata(scCNA)$suggestedK_df
  sug_k <- S4Vectors::metadata(scCNA)$suggestedK

  df <- dplyr::mutate(df, k = as.character(k))

  # df expanded for geom tile
  df_exp <- tidyr::complete(df,
                            k,
                            subclones)

  # common layers
  common_layers <- list(scale_y_discrete(limits = gtools::mixedsort(unique(df$subclones))),
                        scale_x_discrete(limits = gtools::mixedsort(unique(df$k))),
                        theme_classic(), labs(fill = 'jaccard\nsimilarity'))

  if( geom == 'dotplot') {

    p <- ggplot(df, aes(k, subclones)) +
      geom_point(aes(size = bootmean,
                     fill = bootmean),
                 shape = 21) +
      scale_fill_viridis_c(option = 2) +
      common_layers

  }

  if (geom == 'tile') {

    p <- ggplot(df_exp, aes(k, subclones)) +
      geom_tile(aes(fill = bootmean), color = 'black') +
      scale_fill_viridis_c(na.value = 'grey', option = 2) +
      common_layers +
      theme(panel.border = element_rect(fill = NA, size = 3))

  }

  if (geom == 'boxplot') {

    mean_per_k <- df %>%
      dplyr::group_by(k) %>%
      dplyr::summarise(mean_jac = mean(bootmean))

    # adding color for chosen k
    df <- dplyr::mutate(df, chosen = ifelse(df$k == as.character(sug_k),
                                        TRUE,
                                        FALSE))

    p <- ggplot() +
      geom_boxplot(data = df, aes(k, bootmean, fill = chosen),
                   alpha = 1) +
      geom_point(data = mean_per_k, aes(k, mean_jac),
                 fill = 'red',
                 shape = 21,
                 size = 3) +
      scale_fill_manual(values = c("TRUE" = 'khaki', "FALSE" = "grey90")) +
      scale_x_discrete(limits = gtools::mixedsort(unique(df$k))) +
      theme_classic() +
      theme(axis.line.y = element_blank(),
            axis.line.x = element_blank(),
            legend.position = 'none') +
      labs(y = 'jaccard similarity')

  }

  if (geom == 'scatterplot') {

    p <- ggplot(df, aes(x = n_cells, y = bootmean)) +
      geom_point(aes(fill = subclones), shape = 21) +
      stat_smooth(method = 'lm', se = FALSE) +
      facet_wrap(vars(as.numeric(k)), scales = 'free_x') +
      scale_fill_manual(values = subclones_pal(),
                        limits = gtools::mixedsort(unique(df$subclones))) +
      theme_classic() +
      labs(x = "number of cells",
           y = 'jaccard similarity')


  }

  print(p)

}
