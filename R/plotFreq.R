plotFreq <- function(scCNA,
                     threshold = 1,
                     assay = "segment_ratios",
                     label = NULL,
                     geom = c('area', 'line')) {

  geom <- match.arg(geom)

  # gather data
  dat <- as.data.frame(t(SummarizedExperiment::assay(scCNA, assay)))
  meta <- as.data.frame(SummarizedExperiment::colData(scCNA))

  # creating event matrix
  dat_class <-
    as.data.frame(apply(
      dat,
      2,
      cut,
      breaks = c(-Inf, threshold, Inf),
      labels = c("loss", "gain")
    ))

  # if label is provided the dataframe will be split according to the label
  # otherwise use the full dataset
  if (!is.null(label)) {
    meta_vector <- dplyr::pull(meta, label)
    dat_split <- split(dat_class, meta_vector)
  } else {
    dat_split <- dat_class
  }

  # calculating frequency table
  freq_table <- lapply(dat_split, function(x) {

     x %>%
      tidyr::gather(key = 'bin',
                    value = 'value') %>%
      dplyr::mutate(bin = stringr::str_remove(bin,'V')) %>%
      dplyr::mutate(bin = as.numeric(bin)) %>%
      dplyr::group_by(bin) %>%
      dplyr::count(value) %>%
      dplyr::mutate(freq = n/sum(n)) %>%
      dplyr::ungroup() %>%
      tidyr::complete(bin, value, fill = list(freq = 0, n = 0))

  })

  freq_df <- dplyr::bind_rows(freq_table, .id = 'label')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # plot
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (geom == 'area') {
  p <-   ggplot() +
      geom_area(data = subset(freq_df, value == 'gain'),
                aes(bin, freq),
                fill = 'firebrick3') +
      geom_area(data = subset(freq_df, value == 'loss'),
                aes(bin,-freq),
                fill = 'dodgerblue3') +
      facet_wrap(vars(label), ncol = 1)

  }

  if (geom == 'line' && !is.null(label)) {

    p <-  ggplot() +
      geom_line(data = subset(freq_df, value == 'gain'),
                aes(bin, freq, color = label)) +
      geom_line(data = subset(freq_df, value == 'loss'),
                aes(bin,-freq, color = label))

  }

  if (geom == 'line' && is.null(label)) {

    p <-  ggplot() +
      geom_line(data = subset(freq_df, value == 'gain'),
                aes(bin, freq), color = 'firebrick3') +
      geom_line(data = subset(freq_df, value == 'loss'),
                aes(bin,-freq), color = 'dodgerblue3')

  }


}
