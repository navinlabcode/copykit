#' plotAlluvial()
#'
#' Produces an alluvial plot from character elements of the metadata
#'
#' @param scCNA The CopyKit object.
#' @param label A string with two or more  elements from \code{\link[SummarizedExperiment]{colData}}.
#' @param label_colors An optional named vector with the colors of each element
#' from label.
#' @param min_cells An optional numeric to filter stratum that do not reach
#' the minimum amount of cells.
#'
#' @return A ggplot object containing an alluvial plot from ggalluvial
#'
#' @import ggalluvial
#' @import ggplot2
#' @importFrom ggalluvial stat_stratum geom_stratum
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr all_of across group_by count filter
#' @importFrom scales hue_pal
#'
#' @export
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- findClusters(copykit_obj)
#' colData(copykit_obj)$section <- stringr::str_extract(
#'     colData(copykit_obj)$sample,
#'     "(L[0-9]+L[0-9]+|L[0-9]+)"
#' )
#' plotAlluvial(copykit_obj, label = c("subclones", "section"))
plotAlluvial <- function(scCNA,
                         label,
                         label_colors = NULL,
                         min_cells = NULL) {

    # bindings for NSE
    group <- cohort <- NULL

    # thanks for error solving from SO user twedl:
    # https://stackoverflow.com/a/53798038
    StatStratum <- ggalluvial::StatStratum

    meta <- as.data.frame(colData(scCNA))

    # check
    if (all(is.numeric(meta[label]))) {
        stop("label argument must not contain numeric columns.")
    }

    # calculating frequencies across labels
    alluvial_dat <- meta %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(label))) %>%
        dplyr::count() %>%
        ggalluvial::to_lodes_form(
            key = "class",
            value = "group",
            id = "cohort",
            axes = 1:length(label)
        )

    if (!is.null(min_cells)) {
        alluvial_dat <- alluvial_dat %>%
            dplyr::filter(n > min_cells)
    }

    # managing colors
    if (is.null(label_colors)) {
        # defaults
        label_colors <- c(
            superclones_pal(),
            subclones_pal(),
            c(
                "removed" = "#DA614D",
                "kept" = "#5F917A"
            ),
            c(
                "TRUE" = "#396DB3",
                "FALSE" = "#11181D"
            )
        )

        # non defaults
        non_default <- label[label %!in% c(
            "superclones",
            "subclones",
            "is_aneuploid",
            "outlier"
        )]

        non_default_colors <- vector(mode = "list")

        for (i in seq_along(non_default)) {
            # luminescence and brightness
            l <- 65
            h <- 15

            groups_label <- unique(meta[[non_default[i]]])

            non_default_colors[[i]] <-
                structure(scales::hue_pal(
                    h = c(0, 360) + h,
                    l = l
                )(length(groups_label)),
                names = groups_label
                )
            names(non_default_colors)[i] <- non_default[i]

            l <- l - 10
            h <- h + 15
        }

        label_colors <-
            c(label_colors, unlist(unname(non_default_colors)))
    }

    # plot
    p <- ggplot(
        data = alluvial_dat,
        aes(
            x = class,
            stratum = group,
            alluvium = cohort,
            y = n
        )
    ) +
        ggalluvial::geom_flow(
            aes(fill = group),
            stat = "alluvium",
            color = "black",
            alpha = .7,
            width = 1 / 8
        ) +
        ggalluvial::geom_stratum(aes(fill = group), color = "black", width = 1 / 8) +
        geom_text(stat = StatStratum, aes(label = group)) +
        theme_void() +
        theme(
            legend.position = "none",
            axis.text.x = element_text(color = "black")
        ) +
        scale_fill_manual(
            values = label_colors,
            limits = force
        )

    # return plot
    p
}
