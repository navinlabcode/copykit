#' plotMetrics
#'
#' Plots swarm plots from selected elements of \code{\link[SummarizedExperiment]{colData}}.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param metric A character indicating which elements of \code{colData()}
#' should be plotted.
#' @param label A character indicating which element of the \code{colData()} to
#' color the plots.
#' @param dodge.width A numeric that adds dodge between the label categories.
#' @param ncol A Integer specifying the number of columns to be used for the
#' panels of a multi-facet plot.
#'
#'
#' @return A ggplot object with swarm plots of the selected metrics.
#' @export
#' @import ggplot2
#'
#' @examples
#' copykit_obj <- mock_bincounts()
#' copykit_obj <- runMetrics(copykit_obj)
plotMetrics <- function(scCNA,
    metric,
    label = NULL,
    dodge.width = 0,
    ncol = 2) {

    # bindings for NSE
    value <- NULL

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:16:37 2021
    # checks
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:16:44 2021

    # check if metric was supplied
    if (is.null(metric)) {
        stop("Please provide elements of colData as argument to metric.")
    }

    # check if any metric exists
    if (!any(metric %in% names(SummarizedExperiment::colData(scCNA)))) {
        stop("No value from argument metric can be found in colData(scCNA).")
    }
    # check if label exists
    if (!is.null(label)) {
        message("Coloring by: ", label)
    }

    # check if label is length = 1
    if (length(label) > 1) {
        stop("Label argument must be of length = 1.")
    }

    # check if label is an element of colData
    if (!is.null(label) && !(label %in% names(SummarizedExperiment::colData(scCNA)))) {
        stop("Label ", label, " is not a column of the scCNA object.")
    }

    if (any(!(metric %in% names(SummarizedExperiment::colData(scCNA))))) {
        warning(
            "Metric: ",
            paste(metric[!(metric %in% names(SummarizedExperiment::colData(scCNA)))],
                collapse = ", "
            ),
            ", is not a column of the scCNA object."
        )

        metric <- metric[(metric %in% names(SummarizedExperiment::colData(scCNA)))]
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:22:13 2021
    # retrieving data
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:22:18 2021

    df <- as.data.frame(SummarizedExperiment::colData(scCNA))

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:17:41 2021
    # theme setup
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:17:45 2021

    my_theme <- list(
        ggplot2::theme(
            axis.title.x = ggplot2::element_text(colour = "gray28", size = 16),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(colour = "gray28", size = 16),
            axis.text.y = ggplot2::element_text(size = 15),
            axis.line.x = ggplot2::element_blank(),
            legend.position = "right",
            strip.background = element_blank(),
            strip.text = ggplot2::element_text(size = 15),
            legend.text = ggplot2::element_text(size = 16)
        )
    )

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:43:05 2021
    # getting metrics
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 13:44:37 2021

    df_set <- dplyr::select(df, sample, dplyr::all_of(metric)) %>%
        tidyr::gather(
            key = "metric",
            value = "value",
            -sample
        )

    df_label <- dplyr::select(df, sample, !!label) %>%
        tidyr::gather(
            key = "label_name",
            value = "label_value",
            -sample
        )

    df_merge <- dplyr::left_join(df_set, df_label, by = "sample")

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:22:39 2021
    # base plot
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:22:43 2021

    p <- ggplot2::ggplot(df_merge, aes(metric, value)) +
        ggplot2::facet_wrap(vars(metric),
            scales = "free",
            ncol = ncol,
            strip.position = "left"
        ) +
        ggplot2::theme_classic() +
        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 8)) +
        xlab("") +
        ylab("") +
        theme(strip.placement = "outside") +
        my_theme

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:26:06 2021
    # label coloring logic
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 25 14:26:16 2021

    # obtaining label
    if (!is.null(label)) {
        p <- p +
            ggbeeswarm::geom_quasirandom(aes(fill = label_value),
                shape = 21,
                stroke = 0.05,
                dodge.width = dodge.width,
                size = 2.3
            )

        if (label == "superclones") {
            color_lab <- list(ggplot2::scale_fill_manual(
                values = superclones_pal(),
                limits = force
            ))
        } else if (label == "subclones") {
            color_lab <- list(ggplot2::scale_fill_manual(
                values = subclones_pal(),
                limits = force
            ))
        } else if (is.numeric(df_merge$label_value)) {
            color_lab <- list(ggplot2::scale_fill_viridis_c())
        } else {
            color_lab <- list(ggplot2::scale_fill_viridis_d(
                begin = 0.2,
                end = 0.8
            ))
        }

        p <- p + color_lab +
            labs(fill = label)

        # return plot
        p
    } else {
        # else just print the normal without colors
        p <- p +
            ggbeeswarm::geom_quasirandom()

        # return plot
        p
    }
}
