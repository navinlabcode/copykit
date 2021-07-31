#' plotSuggestedK
#'
#' Uses the information from \code{\link{findSugggestedK}} to plot the values
#' of jaccard similarity from the tested k range on \code{\link{findSugggestedK}}.
#'
#' @param scCNA The scCNA object.
#' @param geom A character with the geom to be used for plotting.
#'
#' @details \code{\link{plotSuggestedK}} access the \code{\link[S4Vectors]{metadata}}
#' element suggestedK_df that is saved to the scDNA object after running
#' \code{\link{findSugggestedK}}. The dataframe is used for plotting either a
#' heatmap, when the argument geom = 'tile', or a dotplot when argument geom =
#' 'dotplot' or a boxplot when geom = 'boxplot'.
#'
#' \itemize{
#' \item{geom = 'boxplot':} With geom boxplot the red dots represent the mean
#' jaccard similarity for each assesed k.
#' }
#'
#' @return A ggplot2 object with the plot of different tested k values and their
#' jaccard similarity for each subclone
#' @export
#'
#' @import ggplot2
#' @importFrom tidyr complete
#' @importFrom S4Vectors metadata
#'
#' @examples
plotSuggestedK <- function(scCNA,
                           geom = c('boxplot', 'tile', 'dotplot')) {

  geom = match.arg(geom)

  df <- S4Vectors::metadata(scCNA)$suggestedK_df

  df <- mutate(df, k = as.character(k))

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

    p <- ggplot() +
      geom_boxplot(data = df, aes(k, bootmean),
                   fill = 'grey80') +
      geom_point(data = mean_per_k, aes(k, mean_jac),
                 fill = 'red',
                 shape = 21,
                 size = 3) +
      scale_x_discrete(limits = gtools::mixedsort(unique(df$k))) +
      theme_classic() +
      theme(axis.line.y = element_blank(),
            axis.line.x = element_blank()) +
      labs(y = 'jaccard similarity')

  }

  print(p)

}
