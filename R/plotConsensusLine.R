#' plotConsensusLine
#'
#' @param scCNA The scCNA object
#'
#' @return An interactive plot where different groups
#' can be selected for easier visualization
#' @export
#'
#' @import shiny
#' @import ggplot2
#' @importFrom miniUI miniPage miniContentPanel gadgetTitleBar
#' @importFrom dplyr filter arrange ungroup group_by select row_number
#' @importFrom shiny checkboxGroupInput plotOutput stopApp fillCol
#' @importFrom tidyr gather
#' @importFrom gtools mixedsort
#'
#' @examples
plotConsensusLine <- function(scCNA) {

  ####################
  ## checks
  ####################
  if (nrow(consensus(scCNA)) == 0) {
    stop("Slot consensus is empty. Run calcConsensus()")
  }

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
    ), scale_fill_identity())

  sec_breaks <- c(0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9)
  sec_labels <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)


  # theme
  ggaes <- list(
    scale_x_continuous(
      breaks = xbreaks,
      labels = gsub("chr", "", chrom_rects$chr),
      position = "top",
      expand = c(0, 0),
      sec.axis = sec_axis(
        ~ .,
        breaks = sec_breaks,
        labels = sec_labels,
        name = "genome position (Gb)"
      )
    ),
    theme_classic(),
    theme(
      axis.text.x = element_text(
        angle = 0,
        vjust = .5,
        size = 15
      ),
      axis.text.y = element_text(size = 15),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        size = 1.3
      ),
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 15)
    )
  )
  ####################
  # obtaining and wrangling data
  ####################


  con <- consensus(scCNA)

  con_l <- con %>%
    dplyr::mutate(abspos = chr_ranges$abspos) %>%
    tidyr::gather(key = "group",
                         value = "segment_ratio",
                  -abspos)

  choice <- gtools::mixedsort(unique(con_l$group))

  ####################
  # shiny logic
  ####################

  # tweaks, a list object to set up multicols for checkboxGroupInput
  # alignment thanks to u/Peter
  #https://stackoverflow.com/questions/29738975/how-to-align-a-group-of-checkboxgroupinput-in-r-shiny
  tweaks <-
    list(tags$head(tags$style(HTML("
                                 .multicol {
                                   height: 150px;
                                   -webkit-column-count: 3; /* Chrome, Safari, Opera */
                                   -moz-column-count: 3;    /* Firefox */
                                   column-count: 5;
                                   -moz-column-fill: auto;
                                   -column-fill: auto;
                                 }
                                 "))
    ))

  ui <- miniPage(gadgetTitleBar("Consensus line plot"),
                 miniContentPanel(tweaks,
                                  fillCol(
                                    tags$div(
                                      align = 'left',
                                      class = 'multicol',
                                      checkboxGroupInput(
                                        "checkbox",
                                        label = c(""),
                                        choices = choice,
                                        selected = choice[1]
                                      )
                                    ),

                                    plotOutput("plot", height = "100%"),

                                    # col width
                                    flex = c(1, 2)

                                  )))

  server <- function(input, output, session) {
    # Render the plot
    output$plot <- renderPlot({
      df_plot <- con_l %>%
        dplyr::filter(group %in% input$checkbox)

      p <- ggplot(df_plot) +
        ggchr_back +
        ggaes +
        geom_line(aes(abspos, segment_ratio, color = group),
                  size = 1.2) +
        labs(x = "",
             y = "consensus segment ratio")

      # coloring by superclones or subclones
      if (attr(con, "consensus_by") == "subclones") {
        p <- p + scale_color_manual(values = subclones_pal())
      }

      if (attr(con, "consensus_by") == "superclones") {
        p <- p + scale_color_manual(values = superclones_pal())
      }

      p

    })
    #
    # Handle the Done button being pressed.
    observeEvent(input$done, {
      stopApp(message("Done."))
    })
  }

  runGadget(ui, server)
}
