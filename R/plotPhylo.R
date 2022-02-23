#' plotPhylo()
#'
#' Plots a phylogenetic tree.
#'
#' @author Junke Wang
#'
#' @param scCNA scCNA object.
#' @param label A string with the element of
#' \code{\link[SummarizedExperiment]{colData}}. to annotate the tips of the
#' tree.
#' @param label_colors A named list with colors for the label annotation.
#'  Must match label length
#' @param consensus A boolean indicating if the consensus phylogenetic tree
#' should be plotted.
#' @param group A string that if provided will plot the tip labels as pie charts
#' with the proportions from the provided element from
#' \code{\link[SummarizedExperiment]{colData}}
#'
#' @return A ggplot object with a phylogenetic tree visualization.
#' @export
#'
#' @import ggtree
#' @importFrom scales hue_pal
#' @importFrom dplyr select
#' @importFrom tidyr gather
#' @importFrom scales hue_pal
#' @importFrom ape Ntip
#'
#' @examples
#' set.seed(1000)
#' copykit_obj <- copykit_example_filtered()[,sample(100)]
#' copykit_obj <- findClusters(copykit_obj)
#' copykit_obj <- runPhylo(copykit_obj)
#' plotPhylo(copykit_obj, label = "subclones")
plotPhylo <- function(scCNA,
                      label = NULL,
                      label_colors = NULL,
                      consensus = FALSE,
                      group = NULL) {
    # consensus and not consensus logic
    if (consensus == FALSE) {
        tryCatch(
            phylo(scCNA),
            error = function(e) {
                stop("No phylogeny detected in scCNA object. Please run runPhylo() first.")
            }
        )

        tree <- phylo(scCNA)
    } else {
        tryCatch(
            consensusPhylo(scCNA),
            error = function(e) {
                stop(
                    "No consensus phylogeny detected in scCNA object.
                    Please run runConsensusPhylo() first."
                )
            }
        )

        tree <- consensusPhylo(scCNA)
    }


    # plotting
    if (is.null(label)) {
        # plotting without labels

        ggtree::ggtree(tree, ladderize = FALSE)
    } else {
        # plotting with labels
        if (length(label) > 1) {
            stop("Only one label can be visualized in the tree.")
        }

        metadata <- SummarizedExperiment::colData(scCNA) %>%
            as.data.frame()

        if (!(label %in% colnames(metadata))) {
            stop("Label ", label, " is not a column of colData(scCNA).")
        }

        metadata_anno_df <-
            metadata %>% dplyr::select(dplyr::all_of(label))

        if (consensus == FALSE) {
            metadata_anno_df <- metadata_anno_df[tree$tip.label, , drop = FALSE]
            size <- 1
        }

        if (consensus == TRUE) {
            # Uses the hidden consensus_by attribute from the calcConsensus function
            # to guarantee the same order
            cons_attr <- attr(consensus(scCNA), "consensus_by")

            if (!is.null(group)) {
                # check if group exists
                if (!is.null(group) && !(group %in% colnames(metadata))) {
                    stop("Group ", group, " is not a column of colData(scCNA).")
                }

                # data for groups
                metadata_gr <- metadata %>%
                    droplevels() %>%
                    dplyr::select(!!cons_attr, !!group) %>%
                    tidyr::gather(
                        key = "consensus_var",
                        value = "group_value", -!!cons_attr
                    )

                names(metadata_gr)[1] <- "cons_attr"

                group_split <- metadata_gr$cons_attr
                metadata_gr_list <- split(metadata_gr, group_split)

                # colors for groups
                elements_groups <- sort(unique(as.character(metadata_gr$group_value)))
                n_groups <- length(elements_groups)
                hex <- scales::hue_pal()(n_groups)

                col <- structure(hex,
                    names = elements_groups
                )

                group_colors <- col

                pies <- lapply(metadata_gr_list, function(x) {
                    group_value <- NULL

                    ggplot(
                        x,
                        aes(
                            y = cons_attr,
                            fill = group_value,
                            x = ""
                        )
                    ) +
                        geom_bar(stat = "identity") +
                        scale_fill_manual(values = group_colors, limits = force) +
                        coord_polar(theta = "y", start = 0) +
                        theme_void() +
                        theme(legend.position = "none")
                })

                names(pies) <- names(metadata_gr_list)
            }

            metadata_anno_df <- metadata_anno_df[label] %>%
                dplyr::distinct()

            rownames(metadata_anno_df) <- metadata_anno_df %>%
                dplyr::pull(!!cons_attr)

            metadata_anno_df <-
                metadata_anno_df[names(consensus(scCNA)), , drop = FALSE]

            size <- 5
        }

        if (is.null(label_colors)) {
            # default colors superclones and subclones
            label_colors <- c(
                list(
                    superclones = superclones_pal(),
                    subclones = subclones_pal(),
                    outlier = c(
                        "TRUE" = "#DA614D",
                        "FALSE" = "#5F917A"
                    ),
                    is_aneuploid = c(
                        "TRUE" = "#396DB3",
                        "FALSE" = "#11181D"
                    )
                )
            )

            if (any(grepl(
                paste(c("superclones", "subclones", "outlier", "is_aneuploid"),
                      collapse = "|"),
                label
            ))) {
                # if label is one of the four above, uses the default specified
                # colors above
                label_colors <- label_colors[[label]]
            } else {
                elements <- metadata_anno_df %>%
                    dplyr::pull(label)

                n <- length(elements)
                hues <- seq(15, 375, length = n + 1)
                hex <- scales::hue_pal()(n)

                col <- structure(hex,
                    names = elements
                )

                label_colors <- col
            }
        } else if (is.list(label_colors)) {
            label_colors <- label_colors[[label]]
        }


        list_samples <-
            split(rownames(metadata_anno_df), metadata_anno_df[, label])
        tree <- ggtree::groupOTU(tree, list_samples)
        p <- ggtree::ggtree(tree, ladderize = FALSE) +
            geom_tippoint(aes(color = group), size = size) +
            scale_color_manual(values = label_colors, name = label, limits = force) +
            theme(legend.position = "none") +
            geom_treescale(x = 10)

        if (consensus) {
            p <- p + geom_tiplab(aes(color = group), size = 5, hjust = -0.5)

            if (!is.null(group)) {
                # order of insets must be named by node
                tree_tips <- tree$tip.label[tree$edge[, 2][tree$edge[, 2] <= length(tree$tip.label)]]
                pies <- pies[match(tree_tips, names(pies))]
                names(pies) <- tree$edge[, 2][tree$edge[, 2] <= length(tree$tip.label)]

                p <- ggtree::inset(p,
                    pies,
                    width = 0.06,
                    height = 0.06,
                    hjust = 0
                )
            }
        }

        # return plot
        p
    }
}
