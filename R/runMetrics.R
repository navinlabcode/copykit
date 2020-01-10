#' Run metrics
#'
#' Provides some metrics to quickly assess the quality and overall look of the copy number data.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param n_threads Number of threads used to calculate the distance matrix. Passed to `parallel::mclapply`. As default it uses 1/4 of the detected cores available.
#'
#' @return Adds a the metrics to the scCNA metadata. Those metrics can be used for subsetting the data if desired
#' @return Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @export
#'
#' @examples
#'
runMetrics <- function(scCNA,
                       n_threads = parallel::detectCores() / 4) {

  # checks

  if (!is.numeric(n_threads)) {
    stop("n_threads argument must be numeric.")
  }

  # checks
  if (n_threads > parallel::detectCores()) {
    stop(paste("n_threads argument must be smaller than ", parallel::detectCores()))
  }

  if (n_threads < 1) {
    n_threads = 1
  }


  # theme setup
  my_theme <- list(
    ggplot2::theme(
      axis.title.x = element_text(colour = "gray28", size = 20),
      axis.text.x = element_blank(),
      axis.ticks.x =  element_blank(),
      axis.title.y = element_text(colour = "gray28", size = 20),
      axis.text.y = element_text(size = 15),
      axis.line.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 16)
    )
  )

  ###################
  # Retrieving data

  dat_seg <- copykit::segment_ratios(scCNA)
  dat_rat <- copykit::ratios(scCNA)

  ####################
  # RMSE
  message("Calculating RMSE")
  message(paste("Using", n_threads, "cores."))

  rmses <- parallel::mclapply(1:ncol(dat_seg), function(i) {
    actual <- dat_seg[, i]
    predicted <- dat_rat[, i]
    Metrics::rmse(actual, predicted)
  },
  mc.cores = n_threads)

  rmses <- unlist(rmses)

  SummarizedExperiment::colData(scCNA)$rmse <- rmses

  # rmse plot
  df <- as.data.frame(SummarizedExperiment::colData(scCNA))

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

  if (!is.null(df$pcr_duplicates)) {
    print(cowplot::plot_grid(p1, p2, p3, nrow = 1))
  } else print(p1)

  message("Done.")

}
