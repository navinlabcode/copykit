#' Plot heatmap
#'
#' Plots a heatmap of the copy number data. Each row is a cell and colums represent genomic positions.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param n_threads Number of threads used to calculate the distance. Passed to `amap::Dist`
#'
#' @return A heatmap visualization.
#'
#' @export
#'
#' @examples
#'

plotHeatmap <- function(scCNA,
                        n_threads = 1) {
  #obtaining data
  seg_data <- t(segment_ratios(scCNA))

  # subsetting filtered cells
  if (!is.null(SummarizedExperiment::colData(scCNA)$filtered)) {
    message("Removing filtered out cells.")
    seg_data <- seg_data[SummarizedExperiment::colData(breast_tumor)$filtered == "kept"]
  }

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
      df = as.character(chr[1:nrow(chr), ]),
      show_legend = FALSE,
      show_annotation_name = FALSE,
      which = "column",
      col = list(df = c("1" = "grey88", "2" = "black"))
    )

  if (nrow(seg_data) > 500) {
    message(
      paste(
        "Your dataset has:",
        nrow(seg_data),
        "Cells. Plotting heatmap may take a long time with large number of cells. Set number of threads with n_threads for parallel processing if possible to speed up."
      )
    )
  }

  # ordering cells
  tree <- ape::nj(amap::Dist(seg_data, nbproc = n_threads))

  tree <- ape::ladderize(tree)

  # getting order
  is_tip <- tree$edge[, 2] <= length(tree$tip.label)
  ordered_tips_index <- tree$edge[is_tip, 2]
  tree_tips_order <- tree$tip.label[ordered_tips_index] %>% rev()

  # plotting
  if (is.null(SummarizedExperiment::colData(scCNA)$minor_clusters)) {
    # plotting without clusters

    message("No cluster information detected, use findClusters() to create it.")
    message("Plotting Heatmap.")

    ComplexHeatmap::Heatmap(
      log2(seg_data + 1e-3),
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
    metadata <- colData(scCNA) %>%
      as.data.frame()

    metadata <- metadata[tree_tips_order,]

    metadata_anno_df <- metadata %>%
      dplyr::select(major_clusters,
                    minor_clusters)

    cluster_anno <- ComplexHeatmap::rowAnnotation(df = metadata_anno_df,
                                                  col = list(major_clusters = major_palette,
                                                             minor_clusters = minor_palette),
                                                  show_annotation_name = FALSE)

    seg_data_ordered <- seg_data[tree_tips_order,]

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
