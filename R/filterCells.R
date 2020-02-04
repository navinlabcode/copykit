#' Filter noise cells
#'
#' filterCells uses a k-nearest-neighbor approach to remove cells
#' with random CNA profiles, largely due to noise data.
#' It calculates a correlation matrix and sets a resolution below which non neighbors
#'  will be classified as noise cells.
#'
#' @author Hua-Jun Wu
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param k K-nearest-neighbor, defaults to 5.
#' @param resolution Set's how strict the correlation cut off will be. Defaults to 0.8.
#' @param plot_heatmap Plot a heatmap with annotations of the removed and kept cells (Defaults to TRUE).
#' @param n_threads Number of parallel threads to calculate distances with \code{amap::Dist()}. Defaults to 1/4 of the cores available in your system.
#'
#' @return Adds a filtered cells label to the scCNA metadata. Cells that pass the filtering criteria receive the label "kept", whereas cells that do not pass the filtering criteria receive the label "removed".
#' @return Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @export
#'
#' @examples
#'
#'

filterCells <- function(scCNA,
                        k = 5,
                        resolution = 0.8,
                        plot_heatmap = TRUE,
                        n_threads = parallel::detectCores() / 4) {
  # checks
  if (n_threads < 1) {
    n_threads <- 1
  }

  if (!is.numeric(resolution)) {
    stop("Resolution needs to be a number between 0 and 1")
  }

  if (resolution < 0 || resolution > 1) {
    stop("Resolution needs to be a number between 0 and 1")
  }

  seg <- copykit::segment_ratios(scCNA)

  message("Calculating correlation matrix.")
  dst <- cor(seg)
  dst_knn_df <- apply(as.matrix(dst), 1, function(x) {
    mean(sort(x, decreasing = T)[2:(k + 1)])
  }) %>%
    tibble::enframe(name = "sample",
                    value = "cor")

  dst_knn_df <- dst_knn_df %>%
    dplyr::mutate(filtered = dplyr::case_when(cor >= resolution ~ "kept",
                                              cor < resolution ~ "removed"))

  message(
    "Adding information to metadata. Access with SummarizedExperiment::colData(scCNA)."
  )
  if (identical(SummarizedExperiment::colData(scCNA)$sample,
                dst_knn_df$sample)) {
    SummarizedExperiment::colData(scCNA)$filter_corr_value <- dst_knn_df$cor
    SummarizedExperiment::colData(scCNA)$filtered <- dst_knn_df$filtered

  } else
    stop("Sample names do not match metadata sample info. Check colData(scCNA).")

  if (plot_heatmap == TRUE) {
    message("Plotting heatmap.")

    if (nrow(dst_knn_df) > 500) {
      message(
        paste(
          "Your dataset has:",
          nrow(dst_knn_df),
          "Cells. Plotting heatmap may take a long time with large number of cells.\n Set number of threads with n_threads for parallel processing if possible to speed up."
        )
      )

      message(paste("Using", n_threads, "cores."))

    }

    # filtered annotation
    filter_anno <-
      ComplexHeatmap::rowAnnotation(filtered = dplyr::pull(dst_knn_df, filtered),
                                    col = list(filtered = c(
                                      "kept" = "green2",
                                      "removed" = "firebrick3"
                                    )))

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
        chr_text = ComplexHeatmap::anno_text(v[1:nrow(seg)],
                                             gp = grid::gpar(fontsize = 14)),
        df = as.character(chr[1:nrow(chr), ]),
        show_legend = FALSE,
        show_annotation_name = FALSE,
        which = "column",
        col = list(df = c("1" = "grey88", "2" = "black"))
      )

    # Heatmap
    ht <- ComplexHeatmap::Heatmap(
      log2(t(seg)),
      cluster_rows = function(x) {
        fastcluster::hclust(amap::Dist(x, method = "manhattan", nbproc = n_threads),
                            method = "ward.D2")
      },
      cluster_columns = FALSE,
      use_raster = TRUE,
      raster_device = "CairoPNG",
      top_annotation = chr_bar,
      border = TRUE,
      show_row_names = FALSE,
      show_column_names = FALSE,
      row_split = dplyr::pull(dst_knn_df, filtered),
      left_annotation = filter_anno,
      heatmap_legend_param = list(title = "log2(segratio)")
    )


    print(ht)
  }


  message("Done.")
  return(scCNA)

}
