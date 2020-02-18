#' @title
#' Internal CopyKit functions
#'
#' @description
#' Methods to get or set internal fields from the scCNA class

#' @export
setMethod("segment_ratios", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the segment_ratios data within the assay slot
  SummarizedExperiment::assay(x, "segment_ratios")
})


#' @export
setMethod("ratios", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the ratios data slot
  SummarizedExperiment::assay(x, "ratios")
})

#' @export
setMethod("bin_counts", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the bin_counts data slot
  SummarizedExperiment::assay(x, "bin_counts")
})

#' @export
setMethod("phylo", "scCNA", function(x) {
  # accessor for the phylo slot
      out <- x@phylo
      out
})

#' @export
setReplaceMethod("phylo", "scCNA", function(x, value) {
  # setter method for phylo slot
  x@phylo <- value
  x
})

#' @export
setMethod("distMat", "scCNA", function(x) {
  # accessor for the distMat slot
  out <- x@distMat
  out
})

#' @export
setReplaceMethod("distMat", "scCNA", function(x, value) {
  # setter method for distMat slot
  x@distMat <- value
  x
})

#' @export
setMethod("graph", "scCNA", function(x) {
  # accessor for the getGraph slot
  out <- x@graph
  out
})

#' @export
setReplaceMethod("graph", "scCNA", function(x, value) {
  # setter method for getGraph slot
  x@graph <- value
  x
})

#' @export
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "scCNA", function(object) {
  callNextMethod()
  cat(
    "rowRanges has: ", length(SummarizedExperiment::rowRanges(object)), " ranges\n",
    sep = "",
    "Phylo:"
  )
})


#' @export
# Internal function that plots the interactively heatmap for plotRatioPlot
.interactivelyHeatmap <- function(scCNA) {

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

  #plotting
  ht <- ComplexHeatmap::Heatmap(
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

  print(ht)

  return(seg_data_ordered)

}
