#' plotFreq
#'
#' @param scCNA scCNA object.
#' @param high_threshold A numeric with the threshold above which events are
#'  considered amplifications.
#' @param low_threshold A numeric with the threshold below which events are
#'  considered deletions.
#' @param assay String with the name of the assay to pull data from to plot
#' the frequency plot.
#' @param label A string with the name of the columns from
#' \code{\link[SummarizedExperiment]{colData}} to separate each frequency plot.
#' @param geom A character with the desired geom
#'
#' @details \code{plotFreq} retrieves the data from the desired assay and creates
#' an event matrix based on the high and low thresholds arguments. Values above
#' the high threshold will be classified as gains whereas values below are classified
#' as deletions. The resulting plot is a frequency plot where values above 0
#' represent the frequency of gains and values below 0 represent the frequency of
#' deletions.
#'
#' If the argument 'label' is provided the frequency plot will be calculated
#' separetely for each label. Labels can be any string column from
#' \code{\link[SummarizedExperiment]{colData}}
#'
#' The following geoms are available:
#'
#' \itemize{
#' \item{area}: If geom = 'area' an area plot with the frequency is plotted.
#' If the label argument is provided a different facet will be plotted for each
#' label.
#'
#' \item{line}: If geom = 'line' a line plot with the frequency is plotted.
#' If the label argument lines are overlapped with different colors.
#'
#' }
#'
#' @return A ggplot object with a frequency plot
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr mutate filter group_by ungroup arrange n count bind_rows
#' @importFrom gather complete
#' @importFrom SummarizedExperiment assay rowRanges colData
#'
#' @examples
plotFreq <- function(scCNA,
                     high_threshold = 1.3,
                     low_threshold = 0.7,
                     assay = "segment_ratios",
                     label = NULL,
                     geom = c('area', 'line')) {

  geom <- match.arg(geom)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## aesthetic setup
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # obtaining chromosome lengths
  chr_ranges <-
    as.data.frame(SummarizedExperiment::rowRanges(scCNA))
  chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths


  # obtaining first and last row of each chr
  chr_ranges_start <-  chr_ranges %>%
    dplyr::group_by(seqnames) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(dplyr::row_number() == 1) %>%
    dplyr::ungroup()

  chr_ranges_end <-  chr_ranges %>%
    dplyr::group_by(seqnames) %>%
    dplyr::arrange(seqnames, start) %>%
    dplyr::filter(dplyr::row_number() == dplyr::n()) %>%
    dplyr::ungroup()

  # Creating data frame object for chromosome rectangles shadows
  chrom_rects <- data.frame(
    chr = chr_ranges_start$seqnames,
    xstart = as.numeric(chr_ranges_start$abspos),
    xend = as.numeric(chr_ranges_end$abspos)
  )
  xbreaks <- rowMeans(chrom_rects %>%
                        dplyr::select(xstart,
                                      xend))

  if (nrow(chrom_rects) == 24) {
    chrom_rects$colors <- rep(c("white", "gray"),
                              length(chr_lengths) / 2)
  } else {
    chrom_rects$colors <- c(rep(c("white", "gray"),
                                (length(chr_lengths) / 2)), "white")
  }

  # Creating the geom_rect object
  ggchr_back <-
    list(geom_rect(
      data = chrom_rects,
      aes(
        xmin = xstart,
        xmax = xend,
        ymin = -Inf,
        ymax = Inf,
        fill = colors
      ),
      alpha = .2
    ),
    scale_fill_identity())

  sec_breaks <- c(0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9)
  sec_labels <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)

  # theme
  ggaes <- list(
    scale_x_continuous(
      breaks = xbreaks,
      labels = gsub("chr", "", chrom_rects$chr),
      expand = c(0, 0)
    ),
    theme_classic(),
    xlab("chromosome"),
    ylab("frequency"),
    theme(
      axis.text.x = element_text(
        angle = 0,
        vjust = .5,
        size = 15
      ),
      axis.text.y = element_text(size = 15),
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 15),
      panel.border = element_rect(colour = "black", fill=NA, size = 1.3)
    )
  )

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Data
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  # gather data
  dat <- as.data.frame(t(SummarizedExperiment::assay(scCNA, assay)))
  meta <- as.data.frame(SummarizedExperiment::colData(scCNA))

  # creating event matrix
  dat_class <-
    as.data.frame(apply(
      dat,
      2,
      cut,
      breaks = c(-Inf,low_threshold, high_threshold, Inf),
      labels = c("loss", "neutral", "gain")
    ))

  # if label is provided the dataframe will be split according to the label
  # otherwise use the full dataset
  if (!is.null(label)) {
    meta_vector <- dplyr::pull(meta, label)
    dat_split <- split(dat_class, meta_vector)
  } else {
    dat_split <- list(frequency_plot = dat_class)
  }

  # calculating frequency table
  freq_table <- lapply(dat_split, function(x) {

    colnames(x) <- chr_ranges$abspos

     x %>%
      tidyr::gather(key = 'abspos',
                    value = 'value') %>%
      dplyr::mutate(abspos = as.numeric(abspos)) %>%
      dplyr::group_by(abspos) %>%
      dplyr::count(value) %>%
      dplyr::mutate(freq = n/sum(n)) %>%
      dplyr::ungroup() %>%
      tidyr::complete(abspos, value, fill = list(freq = 0, n = 0))

  })

  freq_df <- dplyr::bind_rows(freq_table, .id = 'label')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # plot
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if (geom == 'area') {
  p <- ggplot() +
    ggchr_back +
    ggaes +
      geom_area(data = subset(freq_df, value == 'gain'),
                aes(abspos, freq),
                fill = 'firebrick3') +
      geom_area(data = subset(freq_df, value == 'loss'),
                aes(abspos,-freq),
                fill = 'dodgerblue3') +
      facet_wrap(vars(label), ncol = 1)

  }

  if (geom == 'line' && !is.null(label)) {

    p <-  ggplot() +
      ggchr_back +
      ggaes +
      geom_line(data = subset(freq_df, value == 'gain'),
                aes(abspos, freq, color = label)) +
      geom_line(data = subset(freq_df, value == 'loss'),
                aes(abspos,-freq, color = label)) +
      theme(legend.position = "bottom")

  }

  if (geom == 'line' && is.null(label)) {

    p <-  ggplot() +
      ggchr_back +
      ggaes +
      geom_line(data = subset(freq_df, value == 'gain') ,
                aes(abspos, freq), color = 'firebrick3') +
      geom_line(data = subset(freq_df, value == 'loss'),
                aes(abspos,-freq), color = 'dodgerblue3')

  }

print(p)

}
