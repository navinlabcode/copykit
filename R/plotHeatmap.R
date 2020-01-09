#' Plot heatmap
#'
#' Plots a heatmap of the copy number data. Each row is a cell and colums represent genomic positions.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @order_cells Methods to order the cells within the heatmap. Accepted values are "phylogeny", "hclust", "graph_search". Defaults to "phylogeny".
#'
#' @return A heatmap visualization.
#'
#' @export
#'
#' @examples
#'

plotHeatmap <- function(scCNA,
                        order_cells = "phylogeny") {
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
    data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], b = chr_rl_c[2:length(chr_rl_c)])
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
      df = as.character(chr[1:nrow(chr),]),
      show_legend = FALSE,
      show_annotation_name = FALSE,
      which = "column",
      col = list(df = c("1" = "grey88", "2" = "black"))
    )

  # ordering cells
  if (order_cells == "phylogeny") {
    tryCatch(
      phylo(scCNA),
      error = function(e) {
        message("No phylogeny detected in scCNA object.")
        runPhylo(scCNA)
      }
    )

    tree <- phylo(scCNA)

    # getting order
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()

    # ordering data
    seg_data_ordered <- seg_data[tree_tips_order, ]

  }

  if (order_cells == "hclust") {
    # checking distance matrix
    if (length(copykit::distMat(scCNA)) == 0) {
      message("No distance matrix detected in the scCNA object.")
      scCNA <-  runDistMat(scCNA, metric = "euclidean")
    }

    if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
      stop(
        "Number of samples in the distance matrix different from number of samples in the scCNA object. Perhaps you filtered your dataset? use copykit::runDistMat() to update it."
      )
    }

    hc <- fastcluster::hclust(distMat(scCNA),
                              method = "ward.D2")

    seg_data_ordered <- seg_data[hc$order, ]

  }

  if (order_cells == "graph_search") {

    if (is.null(SingleCellExperiment::reducedDim(scCNA, 'umap', withDimnames = F))) {
      message("No umap detected, running copykit::runUmap()")
      scCNA <- copykit::runUmap(scCNA)
    }

    tryCatch(
      copykit::graph(scCNA),
      error = function(e) {
        stop("No graph detected. Please run copykit::findClusters()")
      }
    )

    umap_df <-
      SingleCellExperiment::reducedDim(scCNA, 'umap', withDimnames = F)

    g_minor  <- copykit::graph(scCNA)

    g_bs <- igraph::bfs(g_minor, 1)

    g_ord <- as.numeric(g_bs$order)

    seg_data_ordered <- seg_data[g_ord, ]

  }



  # plotting
  if (is.null(SummarizedExperiment::colData(scCNA)$minor_clusters)) {
    # plotting without clusters

    message("No cluster information detected, use findClusters() to create it.")
    message("Plotting Heatmap.")

    ComplexHeatmap::Heatmap(
      log2(seg_data_ordered + 1e-3),
      use_raster = TRUE,
      column_title = "Genomic coordinates",
      column_title_gp = grid::gpar(fontsize = 18),
      column_title_side = "bottom",
      row_title = "Single-Cells",
      row_title_gp = grid::gpar(fontsize = 18),
      heatmap_legend_param = list(title = "log2(segratio)"),
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
    metadata <- SummarizedExperiment::colData(scCNA) %>%
      as.data.frame()

    metadata <- metadata[rownames(seg_data_ordered), ]

    metadata_anno_df <- metadata %>%
      dplyr::select(major_clusters,
                    minor_clusters)

    cluster_anno <-
      ComplexHeatmap::rowAnnotation(
        df = metadata_anno_df,
        col = list(major_clusters = major_palette,
                   minor_clusters = minor_palette),
        show_annotation_name = FALSE
      )

    #plotting
    ComplexHeatmap::Heatmap(
      log2(seg_data_ordered + 1e-3),
      use_raster = TRUE,
      left_annotation = cluster_anno,
      column_title = "Genomic coordinates",
      column_title_gp = grid::gpar(fontsize = 18),
      column_title_side = "bottom",
      row_title = "Single-Cells",
      row_title_gp = grid::gpar(fontsize = 18),
      heatmap_legend_param = list(title = "log2(segratio)"),
      top_annotation = chr_bar,
      cluster_rows = FALSE,
      border = TRUE,
      cluster_columns = FALSE,
      show_column_names = FALSE,
      show_row_names = FALSE,
      show_heatmap_legend = TRUE
    )

  }

}
