#' Plot metrics
#'
#' Provides a plot of the metrics to quickly assess the quality and overall look of the copy number data.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#'
#' @return A plot with the metrics available in the metadata
#' @return Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @export
#' @import ggplot2
#'
#' @examples
#'
#'

plotMetrics <- function(scCNA) {

  # retrieving data
  df <- as.data.frame(SummarizedExperiment::colData(scCNA))

  # check if runMetrics was run
  if (is.null(df$rmse)) {
    message("Metrics not detected.")
    message("Running copykit::runMetrics()")
    scCNA <- runMetrics(scCNA)
  }

  # theme setup
  my_theme <- list(
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(colour = "gray28", size = 20),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x =  ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(colour = "gray28", size = 20),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.line.x = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 16)
    )
  )

  # rmse plot

  p1 <- ggplot2::ggplot(df, aes("RMSE",rmse)) +
    ggbeeswarm::geom_quasirandom() +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    ggplot2::ylab("RMSE") +
    xlab("") +
    my_theme

  # total_reads plot
  if (!is.null(df$total_reads)) {
    p2 <- ggplot2::ggplot(df, aes("tot_reads",total_reads)) +
      ggbeeswarm::geom_quasirandom() +
      ggplot2::theme_classic() +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                                  labels = scales::label_scientific()) +
      ggplot2::ylab("Read Count") +
      xlab("") +
      my_theme
  }

  # pcr_duplicates plot
  if (!is.null(df$pcr_duplicates)) {
    p3 <- ggplot2::ggplot(df, aes("pcr_dups",pcr_duplicates)) +
      ggbeeswarm::geom_quasirandom() +
      ggplot2::theme_classic() +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10),
                                  labels = scales::percent_format()) +
      ggplot2::ylab("PCR Duplicates (%)") +
      xlab("") +
      my_theme
  }

  if (!is.null(df$breakpoint_count)) {
    p4 <- ggplot2::ggplot(df, aes("breakpoint_count",breakpoint_count)) +
      ggbeeswarm::geom_quasirandom() +
      ggplot2::theme_classic() +
      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
      ggplot2::ylab("Breakpoint Count Distribution") +
      xlab("") +
      my_theme
  }

  if (!is.null(df$pcr_duplicates)) {
    print(cowplot::plot_grid(p1, p2, p3, p4, nrow = 2))
  } else print(p1)

}
