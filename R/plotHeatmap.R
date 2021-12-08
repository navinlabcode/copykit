#' plotHeatmap
#'
#' Plots a heatmap of the copy number data.
#' Each row is a cell and colums represent genomic positions.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param assay String with the name of the assay to pull data from to plot heatmap.
#' @param order_cells A string with the desired method to order the cells within
#' @param label A vector with the string names of the columns from
#' \code{\link[SummarizedExperiment]{colData}} for heatmap annotation.
#' @param label_colors A named list with colors for the label annotation.
#'  Must match label length and have the same names as label
#' @param group with the names of the columns from
#' \code{\link[SummarizedExperiment]{colData}} to add a barplot with frequency
#' of the groups to a consensus heatmap.
#' @param row_split A string with the names of the columns from
#' \code{\link[SummarizedExperiment]{colData}} to split the heatmap.
#' @param rounding_error A boolean indicating if the rounding error matrix
#' should be plotted.
#' @param genes A character vector with the HUGO symbol for genes to be annotated
#' on the heatmap.
#' @param consensus A boolean indicating if the consensus heatmap should be
#' plotted.
#' @param n_threads Number of threads passed on to \code{runDistMat}.
#'
#' @return A \code{ComplexHeatmap} object with a heatmap of copy number data
#' where the columns are the genomic positions and each row is a cell.
#'
#' @details
#' \itemize{
#'    \item{order_cells}: If order_cells argument is set to 'consensus_tree'
#'    \code{\link{plotHeatmap}} checks for the existence of a consensus matrix.
#'    From the consensus matrix, a minimum evolution tree is built and cells are
#'     ordered following the order of their respective groups from the tree.
#'     If order_cells is set to 'hclust' cells are ordered according to hierarchical
#'      clustering. 'hclust' calculation can be sped up by changing the parameter
#'      'n_threads' if you have more threads available to use. If order_cells
#'      is set to 'phylogeny' \code{\link{plotHeatmap}} will use the tree stored
#'      in the \code{\link{phylo}} slot of the scCNA object to order the cells.
#'      If the \code{\link{phylo}} is empty \code{\link{plotHeatmap}} will run
#'      \code{\link{runPhylo}} to generate the tree.
#'
#'    \item{label}: A vector with the string names of the columns from
#'    \code{\link[SummarizedExperiment]{colData}} for heatmap annotation. The 'label'
#'    argument can take as many columns as desired as argument as long as they
#'    are elements from \code{\link[SummarizedExperiment]{colData}}.
#'
#'    \item{label_colors}: A named list, list element names must match column
#'    names for \code{\link[SummarizedExperiment]{colData}} and list elements
#'    must match the number of items present in the columns provided in argument
#'     'label'. For example: to set colors for column 'filtered' containing
#'     elements 'kept' or 'filtered' a valid input would be:
#'     'list(filtered = c('kept' = 'green', 'filtered' = 'red))'.
#'      Default colors are provided for 'superclones', 'subclones', 'is_normal',
#'      and 'filtered' that can be overriden with 'label_colors'.
#'
#'    \item{rounding_error}: Must be used with assay = 'integer'.
#'    \code{plotHeatmap} will access the ploidies stored into colData(scCNA)$ploidy
#'     that are generated from \code{\link{calcInteger}} and scale rounded integer
#'      values to the segment means. Later this scaled matrix will be subtracted
#'      from the 'integer' assay from \code{\link{calcInteger}} and the resulting
#'       matrix from this subtraction will be plotted. Useful to visualize regions
#'       of high rounding error. Such regions can indicate issues with the ploidy
#'       scaling in use.
#'
#'    \item{consensus}: If set to TRUE, \code{\link{plotHeatmap}} will search for
#'     the consensus matrix in the slot \code{\link{consensus}} and plot the
#'     resulting matrix. Labels annotations can be added with the argument 'label'.
#'
#' }
#'
#' @seealso \code{\link{calcInteger}} For methods of obtaning the 'integer' assay.
#'
#' @references Zuguang Gu, Roland Eils, Matthias Schlesner, Complex heatmaps
#' reveal patterns and correlations in multidimensional genomic data,
#' Bioinformatics, Volume 32, Issue 18, 15 September 2016, Pages 2847–2849,
#' https://doi.org/10.1093/bioinformatics/btw313
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom circlize colorRamp2
#' @importFrom pals ocean.balance
#' @importFrom S4Vectors metadata
#' @importFrom grDevices colors
#' @importFrom SummarizedExperiment assay
#' @importFrom dplyr select pull all_of mutate group_by
#' @importFrom viridis viridis
#' @importFrom scales hue_pal
#' @importFrom ape Ntip
#' @importFrom stats dist
#' @importFrom tidyr pivot_wider
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' set.seed(1000)
#' copykit_obj <- copykit_obj[, sample(200)]
#' copykit_obj <- findClusters(copykit_obj)
#' copykit_obj <- calcConsensus(copykit_obj)
#' copykit_obj <- runConsensusPhylo(copykit_obj)
#' colData(copykit_obj)$section <- stringr::str_extract(
#'     colData(copykit_obj)$sample,
#'     "(L[0-9]+L[0-9]+|L[0-9]+)"
#' )
#' plotHeatmap(copykit_obj, label = c("section", "subclones"))
plotHeatmap <- function(scCNA,
    assay = "segment_ratios",
    order_cells = c("consensus_tree", "hclust", "phylogeny"),
    label = NULL,
    label_colors = NULL,
    group = NULL,
    consensus = FALSE,
    rounding_error = FALSE,
    genes = NULL,
    row_split = NULL,
    n_threads = 1) {
    # args
    order_cells <- match.arg(order_cells)

    # bindings for NSE
    group_value <- NULL

    # check annotation colors
    if (is.null(label) & !is.null(label_colors)) {
        stop("Please provide a label argument.")
    }

    if (!is.null(label_colors) & !is.list(label_colors)) {
        stop("label_colors argument must be a named list.")
    }

    if (!is.null(label_colors)) {
        if (length(label_colors) != length(label)) {
            stop("Label and Label colors arguments must have the same length.")
        }
    }

    if (!is.null(label) & !is.character(label)) {
        stop("Label must be a character vector.")
    }

    if (is.null(SummarizedExperiment::colData(scCNA)$subclones) &&
        order_cells != "phylogeny") {
        message("Ordering by consensus requires cluster information.")
        message("Switching to hclust.")
        order_cells <- "hclust"
    }

    # rounding error check
    if (rounding_error == TRUE && assay != "integer") {
        stop("Rounding error argument must be used with assay 'integer'.")
    }

    # Uses the hidden consensus_by attribute from the calcConsensus function
    # to plot integer color scheme in case consensus = TRUE
    if (consensus == TRUE) {
        if (attr(consensus(scCNA), "consensus_assay") == "integer") {
            assay <- "integer"
        }
    }


    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 11:57:50 2021
    # setup and colors for assay = integer
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 11:58:01 2021

    if (assay == "integer") {
        # truncate ploidy colors to 2* the mean ploidy
        mean_ploidy <- mean(SummarizedExperiment::colData(scCNA)$ploidy)
        ploidy_trunc <- 2 * round(mean_ploidy)

        # colors
        color_heat <-
            structure(pals::ocean.balance(length(0:ploidy_trunc)),
                names = 0:ploidy_trunc
            )

        # special ploidy colors if ground state rounds to 2
        if (round(mean_ploidy) == 2) {
            color_heat <- structure(
                c(
                    "#3787BA",
                    "#95B8C5",
                    "#F0ECEB",
                    "#D7A290",
                    "#BF583B",
                    "#8D1128",
                    "#3C0912"
                ),
                names = c("0", "1", "2", "3", "4", "5", "6")
            )
        }
    }

    # obtaining data
    seg_data <- t(SummarizedExperiment::assay(scCNA, assay))

    # chromosome bar aesthetic
    chr_ranges <-
        as.data.frame(SummarizedExperiment::rowRanges(scCNA))
    chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths

    if (any(chr_ranges$seqnames == "24") ||
        any(chr_ranges$seqnames == "Y") ||
        any(chr_ranges$seqnames == "chrY")) {
        chr_binary <- rep(c(2, 1), length(chr_lengths) / 2)
    } else {
        chr_binary <- c(rep(c(2, 1), (length(chr_lengths) / 2)), 2)
    }

    chr <-
        data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))

    # getting lengths for chr numbers annotation
    chr_rl_c <- c(1, cumsum(chr_lengths))

    # creating a data frame to calculate rowMeans
    chr_df <-
        data.frame(
            a = chr_rl_c[1:length(chr_rl_c) - 1],
            b = chr_rl_c[2:length(chr_rl_c)]
        )
    chr_l_means <- round(rowMeans(chr_df))

    chrom.names <- c(1:22, "X", "Y")

    # creating the vector for chr number annotations
    v <- vector(length = sum(chr_lengths), mode = "character")
    suppressWarnings(v[chr_l_means] <- chrom.names)
    v[is.na(v)] <- ""

    # chr bar with the chr names
    chr_bar <-
        ComplexHeatmap::HeatmapAnnotation(
            chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data)],
                gp = grid::gpar(fontsize = 14)
            ),
            df = as.character(chr[1:nrow(chr), ]),
            show_legend = FALSE,
            show_annotation_name = FALSE,
            which = "column",
            col = list(df = c("1" = "grey88", "2" = "black"))
        )

    # consensus and not consensus logic
    if (consensus == FALSE) {
        # ordering cells
        if (order_cells == "phylogeny") {
            if (ape::Ntip(phylo(scCNA)) == 0) {
                stop("No phylogeny detected in scCNA object. Use runPhylo")
            }

            if (ape::Ntip(consensusPhylo(scCNA)) == 0) {
                stop("No consensus phylogeny detected in scCNA object.")
            }

            tree <- phylo(scCNA)

            # getting order
            is_tip <- tree$edge[, 2] <= length(tree$tip.label)
            ordered_tips_index <- tree$edge[is_tip, 2]
            tree_tips_order <-
                tree$tip.label[ordered_tips_index] %>% rev()

            # ordering data
            seg_data_ordered <- seg_data[tree_tips_order, ]
        }

        if (order_cells == "hclust") {
            # checking distance matrix
            if (length(copykit::distMat(scCNA)) == 0) {
                message("No distance matrix detected in the scCNA object.")
                scCNA <- runDistMat(scCNA,
                                    metric = "euclidean",
                                    n_threads = n_threads)
            }

            if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
                stop(
                    "Number of samples in the distance matrix different from number
                 of samples in the scCNA object. Perhaps you filtered your
                 dataset?."
                )
            }

            hc <- fastcluster::hclust(distMat(scCNA),
                method = "ward.D2"
            )

            seg_data_ordered <- seg_data[hc$order, ]
        }

        if (order_cells == "consensus_tree") {
            if (nrow(consensus(scCNA)) == 0) {
                scCNA <- calcConsensus(scCNA)
            }

            # metadata info
            consensus_by <- attr(consensus(scCNA), "consensus_by")

            meta <- as.data.frame(colData(scCNA)) %>%
                dplyr::select(sample, !!consensus_by)
            meta_info <- as.character(dplyr::pull(meta, !!consensus_by))

            if (ape::Ntip(consensusPhylo(scCNA)) == 0) {
                stop(
                    "No consensus phylogeny in the CopyKit object."
                )
            } else {
                tree <- consensusPhylo(scCNA)
            }

            # getting order
            is_tip <- tree$edge[, 2] <= length(tree$tip.label)
            ordered_tips_index <- tree$edge[is_tip, 2]
            tree_tips_order <-
                tree$tip.label[ordered_tips_index] %>% rev()

            meta_o <- meta[order(match(meta_info, tree_tips_order)), ]
            seg_data_ordered <- seg_data[meta_o$sample, ]
        }
    } else {
        # getting order from tree
        if (ape::Ntip(consensusPhylo(scCNA)) == 0) {
            stop(
                "Build a consensus with calcConsensus and use runConsensusPhylo
        to store the order of the consensus matrix."
            )
        } else {
            tree <- consensusPhylo(scCNA)
        }

        # getting order
        is_tip <- tree$edge[, 2] <= length(tree$tip.label)
        ordered_tips_index <- tree$edge[is_tip, 2]
        tree_tips_order <-
            tree$tip.label[ordered_tips_index] %>% rev()

        # retrieving data
        seg_data <- as.matrix(t(consensus(scCNA)))

        # ordering data
        seg_data_ordered <- seg_data[tree_tips_order, ]
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:11:34 2021
    # In case assay == integer truncating to a maximum value for ht color
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:11:54 2021

    if (assay == "integer") {
        seg_data_int <- seg_data_ordered
        seg_data_int[seg_data_int > ploidy_trunc] <- ploidy_trunc

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Wed Jun 23 11:55:08 2021
        # Calculating rounding error matrix
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Wed Jun 23 11:55:18 2021

        if (rounding_error == TRUE) {
            # calculate the integer matrix without rounding
            int_nr <-
                as.matrix(SummarizedExperiment::assay(scCNA, "segment_ratios")) %*%
                diag(SummarizedExperiment::colData(scCNA)$ploidy)

            # restoring sample names
            names(int_nr) <-
                names(SummarizedExperiment::assay(scCNA, "segment_ratios"))

            # calculating absolute error and plotting
            err <- abs(int_nr - assay(scCNA, "integer"))

            # making sure order is the same
            err <- err[rownames(seg_data_int)]
        }
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # genes annotation
    #   find_scaffold_genes in internals.R
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (!is.null(genes)) {
        df <- find_scaffold_genes(
            scCNA,
            genes
        )

        mk <-
            ComplexHeatmap::columnAnnotation(
                foo = anno_mark(
                    at = df$pos,
                    labels = df$gene,
                    side = "bottom",
                    labels_gp = grid::gpar(fontsize = 12)
                ),
                show_annotation_name = FALSE,
                show_legend = FALSE
            )
    } else {
        mk <- NULL
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # group bar plot logic for consensus
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (consensus == TRUE) {
        metadata <- as.data.frame(colData(scCNA))

        # Uses the hidden consensus_by attribute from the calcConsensus function
        # to guarantee the same order
        cons_attr <- attr(consensus(scCNA), "consensus_by")

        if (!is.null(group)) {
            # check if group exists
            if (!is.null(group) && !(group %in% colnames(metadata))) {
                stop("Group ", group, " is not a column of colData().")
            }

            # data for groups
            metadata_gr <- metadata %>%
                droplevels() %>%
                dplyr::select(!!cons_attr, !!group) %>%
                tidyr::gather(
                    key = "group_key",
                    value = "group_value", -!!cons_attr
                )

            names(metadata_gr)[1] <- "cons_attr"

            # getting counts for barplot
            metadata_counts <- metadata_gr %>%
                dplyr::group_by(cons_attr) %>%
                dplyr::count(group_value) %>%
                dplyr::mutate(n = n / sum(n)) %>%
                tidyr::pivot_wider(
                    names_from = group_value,
                    values_from = n,
                    id_cols = cons_attr,
                    values_fill = 0
                ) %>%
                as.data.frame()

            rownames(metadata_counts) <- metadata_counts[, 1]
            metadata_counts <- metadata_counts[, -1]

            elements_groups <-
                sort(unique(as.character(metadata_gr$group_value)))

            # ordering
            metadata_counts <- metadata_counts[tree_tips_order, elements_groups]

            # colors for groups
            n_groups <- length(elements_groups)
            hex <- scales::hue_pal()(n_groups)

            col <- structure(hex,
                names = elements_groups
            )

            ha_barplot <-
                rowAnnotation(
                    foo = anno_barplot(metadata_counts,
                        gp = grid::gpar(fill = col)
                    ),
                    show_annotation_name = FALSE
                )
        } else {
            ha_barplot <- NULL
        }
    } else {
        ha_barplot <- NULL
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:12:20 2021
    # Complex heatmap plotting
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:12:32 2021

    message("Plotting Heatmap.")

    # list with arguments for complexheatmap
    complex_args <- list(
        use_raster = TRUE,
        bottom_annotation = mk,
        right_annotation = ha_barplot,
        column_title = "genomic coordinates",
        column_title_gp = grid::gpar(fontsize = 18),
        column_title_side = "bottom",
        row_title = paste0(nrow(seg_data_ordered), " samples"),
        row_title_gp = grid::gpar(fontsize = 18),
        top_annotation = chr_bar,
        cluster_rows = FALSE,
        border = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        show_heatmap_legend = TRUE
    )

    # magick raster cannot plot more than 16k columns.
    # need to remove raster_by_magick on resolutions < 200kb
    # https://githubmemory.com/repo/jokergoo/EnrichedHeatmap/issues/58
    if (ncol(seg_data_ordered) > 16000) {
        complex_args$raster_by_magick <- FALSE
    }

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:12:55 2021
    # Setup for label colData annotation and integer

    # There are two main conditions: if argument label is provided or if
    # the assay contains integer values.
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:13:10 2021

    # plotting
    if (is.null(label)) {
        # plotting without clusters

        if (assay != "integer") {
            suppressMessages(do.call(ComplexHeatmap::Heatmap, c(
                list(
                    matrix = log2(seg_data_ordered + 1e-3),
                    heatmap_legend_param = list(title = "log2 (segratio)")
                ),
                complex_args
            )))
        } else {
            # if assay is integer

            if (rounding_error == FALSE) {
                do.call(ComplexHeatmap::Heatmap, c(
                    list(
                        matrix = seg_data_int,
                        heatmap_legend_param = list(title = "copy number"),
                        col = color_heat
                    ),
                    complex_args
                ))
            } else {
                # if rounding_error == TRUE
                do.call(ComplexHeatmap::Heatmap, c(
                    list(
                        matrix = t(err),
                        heatmap_legend_param = list(title = "rounding error"),
                        col = viridis::viridis(200)
                    ),
                    complex_args
                ))
            }
        }
    } else {
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:05:00 2021
        # If argument label is provided
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Fri Jun 18 12:05:09 2021

        # retrieving metadata
        metadata <- SummarizedExperiment::colData(scCNA) %>%
            as.data.frame()

        if (consensus == FALSE) {
            metadata <- metadata[rownames(seg_data_ordered), ]
        }

        metadata_anno_df <- metadata %>%
            dplyr::select(dplyr::all_of(label))

        if (consensus == TRUE) {
            # Uses the hidden consensus_by attribute from calcConsensus
            # to guarantee the same order
            cons_attr <- attr(consensus(scCNA), "consensus_by")

            # Checks that only the same element of consensus is being
            # used for annotation.
            if (length(label) > 1) {
                stop("Label must be of length 1 for consensus heatmap.")
            }

            if (cons_attr != label) {
                stop(
                    "Consensus heatmap can only be annotated with the same
                     metadata element used for generating the consensus matrix."
                )
            }

            metadata_anno_df <- metadata_anno_df[label] %>%
                dplyr::distinct()

            rownames(metadata_anno_df) <- metadata_anno_df %>%
                dplyr::pull(!!cons_attr)

            metadata_anno_df <-
                metadata_anno_df[rownames(seg_data_ordered), , drop = FALSE]
        }

        if (is.null(label_colors)) {
            # h and l controls the value being picked from the color wheel
            # and the brightness
            h <- 15
            l <- 65

            # cont controls for the continuous scale from viridis package
            cont_options <- c("D", "A", "E", "C", "B")
            cont_i <- 1

            label_colors <- list()

            # default colors superclones and subclones
            label_colors <- c(
                list(
                    superclones = superclones_pal(),
                    subclones = subclones_pal(),
                    filtered = c(
                        "removed" = "#DA614D",
                        "kept" = "#5F917A"
                    ),
                    is_normal = c(
                        "TRUE" = "#396DB3",
                        "FALSE" = "#11181D"
                    )
                ),
                label_colors
            )

            default_labels <-
                c("superclones", "subclones", "filtered", "is_normal")

            for (i in 1:length(label)) {
                if (any(str_detect(
                    label[i],
                    default_labels
                ))) {
                    # if label is one of the four above,
                    # uses the default specified colors above
                    label_colors[i] <-
                        label_colors[default_labels[stringr::str_detect(
                            label[i],
                            default_labels
                        )]]
                    names(label_colors)[i] <- label[i]
                } else if (is.numeric(dplyr::pull(metadata_anno_df, label[i]))) {
                    # if current i metadata element is a numeric vector
                    n <- 300
                    min_v <- min(dplyr::pull(metadata_anno_df, label[i]))
                    max_v <- max(dplyr::pull(metadata_anno_df, label[i]))

                    label_colors[i] <-
                        list(circlize::colorRamp2(
                            seq(min_v, max_v, length = n),
                            viridis::viridis(n, option = cont_options[cont_i])
                        ))
                    names(label_colors)[i] <- label[i]

                    cont_i <- cont_i + 1
                } else {
                    # if current i metadata element is not numeric

                    elements <- metadata_anno_df %>%
                        dplyr::pull(label[i]) %>%
                        unique() %>%
                        as.character() %>%
                        sort()

                    n <- length(elements)
                    hex <- scales::hue_pal(h = c(0, 360) + h, l = 65)(n)

                    col <- structure(hex,
                        names = elements
                    )

                    label_colors[i] <- list(col)
                    names(label_colors)[i] <- label[i]

                    # adding an increment to color on the color wheel and
                    # hue 'brightness' so colors of the next element of the
                    # metadata are not the same as the previous
                    l <- l - 10
                    h <- h + 15
                }
            }
        }

        # removing null elements from the label_colors vector
        # in case they contained elements with default colors
        label_colors[sapply(label_colors, is.null)] <- NULL

        cluster_anno <-
            ComplexHeatmap::rowAnnotation(
                df = metadata_anno_df,
                col = label_colors,
                show_annotation_name = FALSE
            )

        # plotting

        if (!is.null(row_split)) {
            if (length(row_split) > 1) {
                stop("row_split length must be 1")
            } else {
                if (assay != "integer") {
                    do.call(
                        ComplexHeatmap::Heatmap,
                        c(
                            list(
                                matrix = log2(seg_data_ordered + 1e-3),
                                row_split = dplyr::pull(
                                    metadata_anno_df,
                                    row_split
                                ),
                                left_annotation = cluster_anno,
                                heatmap_legend_param = list(title = "log2 (segratio)")
                            ),
                            complex_args
                        )
                    )
                } else {
                    # if assay is integer

                    do.call(
                        ComplexHeatmap::Heatmap,
                        c(
                            list(
                                matrix = seg_data_int,
                                row_split = dplyr::pull(
                                    metadata_anno_df,
                                    row_split
                                ),
                                left_annotation = cluster_anno,
                                heatmap_legend_param = list(title = "copy number"),
                                col = color_heat
                            ),
                            complex_args
                        )
                    )
                }
            }
        } else {
            if (assay != "integer") {
                do.call(ComplexHeatmap::Heatmap, c(
                    list(
                        matrix = log2(seg_data_ordered + 1e-3),
                        left_annotation = cluster_anno,
                        heatmap_legend_param = list(title = "log2 (segratio)")
                    ),
                    complex_args
                ))
            } else {
                # if assay == integer and rounding_error == FALSE it will plot
                # the integer matrix
                # otherwise it will plot the rounding error matrix

                if (rounding_error == FALSE) {
                    do.call(ComplexHeatmap::Heatmap, c(
                        list(
                            matrix = seg_data_int,
                            left_annotation = cluster_anno,
                            heatmap_legend_param = list(title = "copy number"),
                            col = color_heat
                        ),
                        complex_args
                    ))
                } else {
                    # if rounding_error == TRUE
                    do.call(ComplexHeatmap::Heatmap, c(
                        list(
                            matrix = t(err),
                            left_annotation = cluster_anno,
                            heatmap_legend_param = list(title = "rounding error"),
                            col = viridis::viridis(200)
                        ),
                        complex_args
                    ))
                }
            }
        }
    }
}
