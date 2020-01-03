#' Plot ratio plot
#'
#' plotRatioPlot allows for a visualization of the segment ratios together with the ratios.
#' It is helpful to observe the fit of the segmentation to the data.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param sample_name character vector with the name of the sample to be visualized
#'
#' @return A violin plot with the segment ratios for the genes of interest.
#' @import BiocGenerics
#' @import S4Vectors
#'
#' @export
#'
#' @examples
#'

plotRatioPlot <- function(scCNA,
                          sample_name) {


  ####################
  ## aesthetic setup
  ####################

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
  chrom_rects <- data.frame(chr = chr_ranges_start$seqnames,
                            xstart = as.numeric(chr_ranges_start$abspos),
                            xend = as.numeric(chr_ranges_end$abspos))
  xbreaks <- rowMeans(chrom_rects %>%
                        dplyr::select(xstart,
                                      xend))

  if (length(chr_rects == 24)) {
    chrom_rects$colors <- rep(c("white", "gray"), length(chr_lengths) / 2)
  } else {
    chrom_rects$colors <- c(rep(c("white", "gray"), (length(chr_lengths) / 2)), "white")
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
    ), scale_fill_identity())

  sec_breaks <- c(0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9)
  sec_labels <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)

  # theme
  ggaes <- list(
    scale_x_continuous(breaks = xbreaks,
                       labels = gsub("chr", "", chrom_rects$chr),
                       position = "top",
                       expand = c(0,0),
                       sec.axis = sec_axis(~., breaks = sec_breaks, labels = sec_labels, name = "Genome Position (Gb)")),
    theme_classic(),
    theme(
      axis.text.x = element_text(angle = 0,
                                 vjust = .5,
                                 size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.y.right = element_text(margin = margin(l = 10)),
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 15)
    )
  )

  ###############
  ## Data setup
  ###############

  abspos <- chr_ranges$abspos

  dat_seg <- copykit::segment_ratios(scCNA) %>%
    dplyr::mutate(abspos = abspos)

  dat_ratios <- copykit::ratios(scCNA) %>%
    dplyr::mutate(abspos = abspos)


  # ggplot will need the long format tables, using gather
  dat_ratios_l <-
    tidyr::gather(
      data = dat_ratios,
      key = "sample",
      value = "ratio",
      -abspos
    )
  dat_seg_l <-
    tidyr::gather(
      data = dat_seg,
      key = "sample",
      value = "segment_ratio",
      -abspos
    )

  if (nrow(dat_ratios_l) == nrow(dat_seg_l)) {
    df <- dat_ratios_l %>%
      dplyr::mutate(segment_ratio = dat_seg_l$segment_ratio)
  } else stop("Nrow in copykit::segment_ratios() assay different than nrow in copykit::ratios().")

  ###############
  ## Plot
  ###############

  p <- ggplot(df %>% filter(sample == sample_name)) +
    ggchr_back +
    ggaes +
    geom_point(aes(abspos, log2(ratio+1e-3)),
               shape = 20,
               col = "gray",
               size = 1,
               alpha = .7) +
    geom_line(aes(abspos, log2(segment_ratio+1e-3)), col = "black",
              size = 1.2) +
    xlab("") +
    ylab("Log2(Ratios)") +
    ggtitle(paste(toupper(sample_name)))

  print(p)

}
