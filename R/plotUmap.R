#' plotUmap
#'
#' Plots UMAP embedding stored in \code{\link[SingleCellExperiment]{reducedDim}}
#' slot.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA The CopyKit object.
#' @param embedding String with the name of the reducedDim to pull data from.
#' @param label A string with the elements from
#' \code{\link[SummarizedExperiment]{colData}} to color the umap points.
#'
#' @details A reduced dimension representation with UMAP in the slot
#' \code{\link[SingleCellExperiment]{reducedDim}} from the scCNA object.
#'
#' Columns from \code{\link[SummarizedExperiment]{colData}} can
#' be used as an argument for 'label' to color the points on the plot.
#'
#' @return A ggplot object containing the reduced dimensions UMAP plot.
#'
#' @export
#'
#' @importFrom ggnewscale new_scale_color
#' @import ggplot2
#' @examples
#' set.seed(1000)
#' copykit_obj <- copykit_example_filtered()[,sample(300)]
#' copykit_obj <- runUmap(copykit_obj)
#'
#' plotUmap(copykit_obj)
#'
#' copykit_obj <- findClusters(copykit_obj)
#'
#' plotUmap(copykit_obj, label = "subclones")
#'
#' colData(copykit_obj)$section <- stringr::str_extract(
#'     colData(copykit_obj)$sample,
#'     "(L[0-9]+L[0-9]+|L[0-9]+)"
#' )
#'
#' plotUmap(copykit_obj, label = c("section"))
plotUmap <- function(scCNA,
                     embedding = "umap",
                     label = NULL) {

    # bindings for NSE objects
    V1 <- V2 <- NULL

    message("Plotting Umap.")

    # retrieving data
    df <- as.data.frame(SummarizedExperiment::colData(scCNA))
    umap_df <- SingleCellExperiment::reducedDim(scCNA, embedding) %>%
        as.data.frame()

    # check if label exists
    if (!is.null(label) && !(label %in% colnames(df))) {
        stop("Label ", label, " is not a column of the scCNA object.")
    }

    if (!is.null(label)) {
        message("Coloring by: ", label, ". ")
    }

    # theme setup
    my_theme <- list(
        ggplot2::theme(
            axis.title.x = element_text(size = 14),
            axis.text.x = element_text(size = 12),
            axis.title.y = element_text(size = 14),
            axis.text.y = element_text(size = 12),
            axis.line = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            legend.position = "right",
            legend.text = element_text(size = 14)
        ),
        xlab("umap1"),
        ylab("umap2")
    )

    # Base plot
    p <- ggplot(umap_df, aes(V1, V2)) +
        theme_classic() +
        my_theme

    if (is.null(label)) {
        p <- p +
            geom_point()

        # return plot
        return(p)
    }

    if (all(label == "subclones")) {
        p <- p +
            geom_point(aes(fill = as.factor(
                SummarizedExperiment::colData(scCNA)$subclones
            )),
            size = 2.5,
            shape = 21,
            stroke = 0.1
            ) +
            scale_fill_manual(
                values = subclones_pal(),
                name = "subclones",
                limits = force
            )

        # return plot
        return(p)
    }

    if (all(label == "superclones")) {
        p <- p +
            geom_point(aes(fill = as.factor(
                SummarizedExperiment::colData(scCNA)$superclones
            )),
            size = 2.5,
            shape = 21
            ) +
            scale_fill_manual(
                values = superclones_pal(),
                name = "superclones",
                limits = force
            )

        # return plot
        return(p)
    }

    if ("subclones" %in% label && "superclones" %in% label) {
        p <- p +
            geom_point(
                aes(
                    x = V1,
                    y = V2,
                    color = SummarizedExperiment::colData(scCNA)$superclones
                ),
                alpha = 1,
                size = 5
            ) +
            scale_color_manual(
                values = superclones_pal(),
                name = "superclones",
                limits = force
            ) +
            ggnewscale::new_scale_color() +
            geom_point(aes(
                x = V1,
                y = V2,
                fill = as.factor(SummarizedExperiment::colData(scCNA)$subclones)
            ),
            size = 2.5,
            shape = 21,
            stroke = 0.1
            ) +
            scale_fill_manual(
                values = subclones_pal(),
                name = "subclones",
                limits = force
            )

        # return plot
        return(p)
    }

    if (!is.null(label) && !("subclones" %in% label && "superclones" %in% label)) {
        if (length(label) > 1) {
            stop("Label must be of length 1.")
        }

        lab <- dplyr::pull(df,
            var = label
        )

        p <- p +
            geom_point(aes(fill = lab),
                size = 2.5,
                shape = 21,
                stroke = 0.1
            ) +
            theme_classic() +
            labs(fill = label) +
            my_theme

        # coloring by continuos variable
        if (is.numeric(lab)) {
            p <- p +
                geom_point(aes(fill = lab),
                    size = 2.5,
                    shape = 21,
                    stroke = 0.1
                ) +
                ggplot2::scale_fill_viridis_c()
        }

        # return plot
        p
    }
}
