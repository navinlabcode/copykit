#' Plot heatmap
#'
#' Plots a heatmap of the copy number data. Each row is a cell and colums represent genomic positions.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param order_cells Methods to order the cells within the heatmap. Accepted values are "phylogeny", "hclust", "graph_search". Defaults to "hclust".
#' @param label Character. Annotate heatmap by an element of metadata. Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @param label_colors List. Named list with colors for the label annotation. Must match label length
#' @param row_split Character. Element of the metadata to split the heatmap. Must have length = 1.
#' @param consensus Boolean. Indicates if the consensus heatmap should be plotted.
#' @param use_default_colors Boolean. Use the default colors when annotating superclones, subclones or filtering.
#'
#' @return A heatmap visualization.
#'
#' @export
#'
#' @import ComplexHeatmap
#' @importFrom dplyr select pull
#' @examples
#'

plotHeatmap <- function(scCNA,
                        order_cells = "consensus_tree",
                        label = NULL,
                        label_colors = NULL,
                        use_default_colors = TRUE,
                        consensus = FALSE,
                        row_split = NULL) {
  # check annotation colors
  if (is.null(label) & !is.null(label_colors)) {
    stop("Please provide a label argument if colors are being specified for it.")
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

  #obtaining data
  seg_data <- t(segment_ratios(scCNA))

  #chromosome bar aesthetic
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
    data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1],
               b = chr_rl_c[2:length(chr_rl_c)])
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
                                           gp = grid::gpar(fontsize = 14)),
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
      tryCatch(
        phylo(scCNA),
        error = function(e) {
          message("No phylogeny detected in scCNA object.")
        },
        finally = {
          scCNA <- runPhylo(scCNA)
        }
      )

      tree <- phylo(scCNA)

      # getting order
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()

      # ordering data
      seg_data_ordered <- seg_data[tree_tips_order,]

    }

    if (order_cells == "hclust") {
      # checking distance matrix
      if (length(copykit::distMat(scCNA)) == 0) {
        message("No distance matrix detected in the scCNA object.")
        scCNA <-  runDistMat(scCNA, metric = "euclidean")
      }

      if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
        stop(
          "Number of samples in the distance matrix different
          from number of samples in the scCNA object.
        Perhaps you filtered your dataset? use copykit::runDistMat() to update it."
        )
      }

      hc <- fastcluster::hclust(distMat(scCNA),
                                method = "ward.D2")

      seg_data_ordered <- seg_data[hc$order,]

    }

    if (order_cells == "consensus_tree") {

      if (is.null(SummarizedExperiment::colData(scCNA)$subclones)) {
        stop("Ordering by consensus requires cluster information. use findClusters(scCNA)")
      }

      if (nrow(consensus(scCNA)) == 0) {
        scCNA <- calcConsensus(scCNA)
        scCNA <- runConsensusPhylo(scCNA)
      }

      # metadata info
      consensus_by <- attr(consensus(scCNA), "consensus_by")

      meta <- as.data.frame(colData(scCNA)) %>%
        dplyr::select(sample, !!consensus_by)
      meta_info <- as.character(dplyr::pull(meta, !!consensus_by))

      tree <- consensusPhylo(scCNA)

      # getting order
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()

      meta_o <- meta[order(match(meta_info, tree_tips_order)), ]
      seg_data_ordered <- seg_data[meta_o$sample,]

    }

    if (order_cells == "graph_search") {
      if (is.null(SingleCellExperiment::reducedDim(scCNA,
                                                   'umap',
                                                   withDimnames = F))) {
        message("No umap detected, running copykit::runUmap()")
        scCNA <- copykit::runUmap(scCNA)
      }

      tryCatch(
        copykit::graph(scCNA),
        error = function(e) {
          stop("No graph detected. Please run copykit::findClusters()")
        }
      )

      if (igraph::gorder(graph(scCNA)) != ncol(copykit::segment_ratios(scCNA))) {
        stop(
          "Number of vertices in graph different than number of samples.
           Please run copykit::findClusters() and try again."
        )
      }

      umap_df <-
        SingleCellExperiment::reducedDim(scCNA, 'umap', withDimnames = F)

      g_minor  <- copykit::graph(scCNA)

      g_bs <- igraph::bfs(g_minor, 1)

      g_ord <- as.numeric(g_bs$order)

      seg_data_ordered <- seg_data[g_ord,]

    }

  } else {
    seg_data_ordered <- as.data.frame(t(consensus(scCNA)))
  }


  # plotting
  if (is.null(label)) {
    # plotting without clusters

    message("Plotting Heatmap.")

    ComplexHeatmap::Heatmap(
      log2(seg_data_ordered + 1e-3),
      use_raster = TRUE,
      column_title = "genomic coordinates",
      column_title_gp = grid::gpar(fontsize = 18),
      column_title_side = "bottom",
      row_title = paste0(nrow(seg_data_ordered), " samples"),
      row_title_gp = grid::gpar(fontsize = 18),
      heatmap_legend_param = list(title = "log2 (segratio)"),
      top_annotation = chr_bar,
      cluster_rows = FALSE,
      border = TRUE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = FALSE,
      show_heatmap_legend = TRUE
    )
  } else {
    #cluster annotation

    # retrieving metadata
    metadata <- SummarizedExperiment::colData(scCNA) %>%
      as.data.frame()

    if (consensus == FALSE) {
      metadata <- metadata[rownames(seg_data_ordered), ]
    }

    metadata_anno_df <- metadata %>%
      dplyr::select(dplyr::all_of(label))

    if (consensus == TRUE) {

      # Uses the hidden consensus_by attribute from the calcConsensus function
      # to guarantee the same order
      cons_attr <- attr(consensus(scCNA), "consensus_by")

      metadata_anno_df <- metadata_anno_df[label] %>%
        dplyr::distinct()

      rownames(metadata_anno_df) <- metadata_anno_df %>%
        dplyr::pull(!!cons_attr)

      metadata_anno_df <- metadata_anno_df[names(consensus(scCNA)),,drop = FALSE]

    }

    if (use_default_colors == TRUE) {
      label_colors <- c(
        list(superclones = superclones_pal(),
             subclones = subclones_pal(),
             filtered = c("kept" = "green2",
                          "removed" = "firebrick3"),
             is_normal = c("FALSE" = "#F8766D",
                           "TRUE" = "#00BFC4")),
        label_colors
      )
    }

    if (is.null(label_colors)) {
      cluster_anno <-
        ComplexHeatmap::rowAnnotation(df = metadata_anno_df,
                                      show_annotation_name = FALSE)
    } else {
      cluster_anno <-
        ComplexHeatmap::rowAnnotation(df = metadata_anno_df,
                                      col = label_colors,
                                      show_annotation_name = FALSE)
    }

    #plotting
    complex_args <- list(
      use_raster = TRUE,
      left_annotation = cluster_anno,
      column_title = "genomic coordinates",
      column_title_gp = grid::gpar(fontsize = 18),
      column_title_side = "bottom",
      row_title = paste0(nrow(seg_data_ordered), " samples"),
      row_title_gp = grid::gpar(fontsize = 18),
      heatmap_legend_param = list(title = "log2 (segratio)"),
      top_annotation = chr_bar,
      cluster_rows = FALSE,
      border = TRUE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = FALSE,
      show_heatmap_legend = TRUE
    )


    if (!is.null(row_split)) {
      if (length(row_split) > 1) {
        stop("row_split length must be 1")
      } else {
        do.call(ComplexHeatmap::Heatmap,
                c(
                  list(
                    matrix = log2(seg_data_ordered + 1e-3),
                    row_split = dplyr::pull(metadata_anno_df, row_split)
                  ),
                  complex_args
                ))
      }
    } else {
      do.call(ComplexHeatmap::Heatmap, c(list(matrix = log2(
        seg_data_ordered + 1e-3
      )), complex_args))
    }

  }

}
