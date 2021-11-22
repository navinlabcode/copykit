#' @title
#' Internal CopyKit functions
#'
#' @description
#' Methods to get or set internal fields from the CopyKit class

#' @export
setMethod("segment_ratios", "CopyKit", function(x, withDimnames = TRUE) {
  # accessor for the segment_ratios data within the assay slot
  SummarizedExperiment::assay(x, "segment_ratios")
})


#' @export
setMethod("ratios", "CopyKit", function(x, withDimnames = TRUE) {
  # accessor for the ratios data slot
  SummarizedExperiment::assay(x, "ratios")
})

#' @export
setMethod("bincounts", "CopyKit", function(x, withDimnames = TRUE) {
  # accessor for the bincounts data slot
  SummarizedExperiment::assay(x, "bincounts")
})

#' @export
setMethod("consensus", "CopyKit", function(x, withDimnames = TRUE) {
  # accessor for the consensus data slot
  out <- x@consensus
  out
})

#' @export
setReplaceMethod("consensus", "CopyKit", function(x, value) {
  # setter method for phylo slot
  x@consensus <- value
  x
})

#' @export
setMethod("phylo", "CopyKit", function(x) {
  # accessor for the phylo slot
  out <- x@phylo
  out
})

#' @export
setReplaceMethod("phylo", "CopyKit", function(x, value) {
  # setter method for phylo slot
  x@phylo <- value
  x
})

#' @export
setMethod("consensusPhylo", "CopyKit", function(x) {
  # accessor for the consensusPhylo slot
  out <- x@consensusPhylo
  out
})

#' @export
setReplaceMethod("consensusPhylo", "CopyKit", function(x, value) {
  # setter method for consensusPhylo slot
  x@consensusPhylo <- value
  x
})

#' @export
setMethod("distMat", "CopyKit", function(x) {
  # accessor for the distMat slot
  out <- x@distMat
  out
})

#' @export
setReplaceMethod("distMat", "CopyKit", function(x, value) {
  # setter method for distMat slot
  x@distMat <- value
  x
})

#' @export
setMethod("graph", "CopyKit", function(x) {
  # accessor for the getGraph slot
  out <- x@graph
  out
})

#' @export
setReplaceMethod("graph", "CopyKit", function(x, value) {
  # setter method for getGraph slot
  x@graph <- value
  x
})

#' @export
#' @importFrom ape Ntip Nnode
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "CopyKit", function(object) {
  callNextMethod()
  cat(
    "rowRanges has: ",
    length(SummarizedExperiment::rowRanges(object)),
    " ranges\n",
    sep = "",
    "Phylo: ",
    "Phylogenetic tree with ",
    ape::Ntip(phylo(object)),
    " tips and ",
    ape::Nnode(phylo(object)),
    " nodes",
    "\n",
    "consensus dim: ",
    nrow(consensus(object)),
    " ",
    ncol(consensus(object))
  )
})

#' @export
`%!in%` <- Negate(`%in%`)

# suppress cat output from function
# thanks to https://stackoverflow.com/a/54136863
#' @export
.quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

# inverse freeman-tukey transformation
# https://www.biorxiv.org/content/10.1101/2020.06.08.140673v2.full.pdf
#' @export
#' @keywords internal
.invft <- function(y) (y^2-1)^2/(4*y^2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# color palettes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' @export
#superclones palette
superclones_pal <- function() {
  structure(
    c(
      "#66C5CC",
      "#F6CF71",
      "#F89C74",
      "#DCB0F2",
      "#87C55F",
      "#9EB9F3",
      "#FE88B1",
      "#C9DB74",
      "#8BE0A4",
      "#B497E7",
      "#D3B484",
      "#B3B3B3",
      "#88CCEE",
      "#CC6677",
      "#DDCC77",
      "#117733",
      "#332288",
      "#AA4499",
      "#44AA99",
      "#999933",
      "#882255",
      "#661100",
      "#6699CC",
      "#888888"
    ),
    names = paste0('s', 1:24)
  )
}

#' @export
#subclones palette
subclones_pal <- function() {
  structure(
    c("gray",
      "#5050FF",
      "#CE3D32",
      "#749B58",
      "#F0E685",
      "#466983",
      "#BA6338",
      "#5DB1DD",
      "#802268",
      "#6BD76B",
      "#D595A7",
      "#924822",
      "#837B8D",
      "#C75127",
      "#D58F5C",
      "#7A65A5",
      "#E4AF69",
      "#3B1B53",
      "#CDDEB7",
      "#612A79",
      "#AE1F63",
      "#E7C76F",
      "#5A655E",
      "#CC9900",
      "#99CC00",
      "#A9A9A9",
      "#CC9900",
      "#99CC00",
      "#33CC00",
      "#00CC33",
      "#00CC99",
      "#0099CC",
      "#0A47FF",
      "#4775FF",
      "#FFC20A",
      "#FFD147",
      "#990033",
      "#991A00",
      "#996600",
      "#809900",
      "#339900",
      "#00991A",
      "#009966",
      "#008099",
      "#003399",
      "#1A0099",
      "#660099",
      "#990080",
      "#D60047",
      "#FF1463",
      "#00D68F",
      "#14FFB1"
    ),
    names = paste0('c', 0:51)
  )
}

#' @author Darlan Conterno Minussi
#' @export
#' @internal
find_scaffold_genes <- function(scCNA,
                                genes) {

  # genome assembly
  if (S4Vectors::metadata(scCNA)$genome == "hg19") {
    genes_assembly <- hg19_genes
  }

  if (S4Vectors::metadata(scCNA)$genome == "hg38") {
    genes_assembly <- hg38_genes
  }

  # getting ranges from scCNA object
  ranges <- SummarizedExperiment::rowRanges(scCNA)

  # subsetting to only the desired genes
  genes_features <- BiocGenerics:::subset(genes_assembly,
                                          symbol %in% genes)

  all_genes <- genes_assembly$symbol %>%
    unlist() %>%
    unname()

  missing_genes <- genes[genes %!in% all_genes]

  if (!rlang::is_empty(missing_genes)) {
    warning(
      base::paste(
        "Genes:",
        paste(missing_genes,
              collapse = ", "),
        ",could not be found. Maybe you need to use a different gene alias?"
      )
    )
  }

  #finding overlaps
  olaps <-
    suppressWarnings(GenomicRanges::findOverlaps(genes_features,
                                                 ranges,
                                                 ignore.strand = TRUE))

  # creating a data frame that will contain the genes and positions (index) in the
  # pipeline ranges.
  # some genes might overlap more than one range (more than one bin), in this case
  # only one will be kept
  df <-
    tibble::tibble(gene = as.character(genes_features$symbol[S4Vectors::queryHits(olaps)]),
                   pos = S4Vectors::subjectHits(olaps)) %>%
    dplyr::distinct(gene, .keep_all = TRUE)

  # checking for genes that might have been excluded from the varbin pipeline
  blk_list <- genes[genes %!in% missing_genes]
  blk_list <- blk_list[blk_list %!in% df$gene]

  if (!rlang::is_empty(blk_list)) {
    warning(
      base::paste(
        "Genes:",
        paste(blk_list,
              collapse = ", "),
        "are in excluded regions of the Varbin pipeline and can't be plotted."
      )
    )
  }

  return(df)

}


#' @author Alexander Davis
#' @export
l2e.normal.sd <- function(xs)
{
  # Need at least two values to get a standard deviation
  stopifnot(length(xs) >= 2)
  optim.result <- stats::optimize(
    # L2E loss function
    f=function(sd)
      # "Data part", the sample average of the likelihood
      -2 * mean(stats::dnorm(xs, sd=sd)) +
      # "Theta part", the integral of the squared density
      1/(2*sqrt(pi)*sd),
    # Parameter: standard deviation of the normal distribution fit
    interval = c(0, diff(range(xs))))
  return(optim.result$minimum)
}

# A function for estimating the index of dispersion, which is used when
# estimating standard errors for each segment mean

#' @author Alexander Davis
#' @export
overdispersion <- function(v)
{
  # 3 elements, 2 differences, can find a standard deviation
  stopifnot(length(v) >= 3)
  # Differences between pairs of values
  y <- v[-1]
  x <- v[-length(v)]
  # Normalize the differences using the sum. The result should be around zero,
  # plus or minus square root of the index of dispersion
  vals.unfiltered <- (y-x)/sqrt(y+x)
  # Remove divide by zero cases, and--considering this is supposed to be count
  # data--divide by almost-zero cases
  vals <- vals.unfiltered[y + x  >= 1]
  # Check that there's anything left
  stopifnot(length(vals) >= 2)
  # Assuming most of the normalized differences follow a normal distribution,
  # estimate the standard deviation
  val.sd <- l2e.normal.sd(vals)
  # Square this standard deviation to obtain an estimate of the index of
  # dispersion
  iod <- val.sd^2
  # subtract one to get the overdispersion criteria
  iod.over <- iod -1
  # normalizing by mean bincounts
  iod.norm <- iod.over/mean(v)
  return(iod.norm)

}

#' @author Junke Wang
#' @export
parCor <- function(x, BPPARAM=BiocParallel::bpparam())
{
  ncol <- ncol(x)

  ## skip parallelization if # of cell less than 2000
  if (ncol < 2001) {
    return(cor(x))
  }

  nSplit <- floor(ncol/1000)
  lSplit <- floor(ncol / nSplit)
  iSplit <- vector("list", nSplit)
  for (i in seq_len((nSplit - 1))) {
    iSplit[[i]] <- (lSplit * (i - 1) + 1):(lSplit * i)
  }
  iSplit[[nSplit]] <- (lSplit * (nSplit - 1) + 1):ncol

  comb <- expand.grid(1:nSplit, 1:nSplit)
  comb <- unique(t(apply(comb, 1, sort)))

  ## result matrix
  result <- BiocParallel::bplapply(X = 1:nrow(comb),
                                   FUN = function(i){
                                     a <- list()
                                     a[[1]] <- iSplit[[comb[i,1]]]
                                     a[[2]] <- iSplit[[comb[i,2]]]
                                     a[[3]] <- cor(x[,a[[1]]], x[,a[[2]]])
                                     a
                                   },
                                   BPPARAM = BPPARAM)

  res_parcor_reserve <- matrix(, nrow = ncol, ncol = ncol)
  for(i in 1:length(result)){
    res_parcor_reserve[result[[i]][[1]],result[[i]][[2]]] <- result[[i]][[3]]
    res_parcor_reserve[result[[i]][[2]],result[[i]][[1]]] <- t(result[[i]][[3]])
  }

  colnames(res_parcor_reserve) <- colnames(x)
  rownames(res_parcor_reserve) <- colnames(x)
  return(res_parcor_reserve)
}

#' @author Darlan Conterno Minussi
#' @export
#' @importFrom  S4Vectors metadata
#' @importFrom rlang chr
#' @internal
copykit_example <- function() {
  #ranges
  hg38_rg_edit <- hg38_rg[, -c(4:5)]
  hg38_rg_editnoY <- dplyr::filter(hg38_rg_edit, chr != "chrY")
  rg_hg38 <- GenomicRanges::makeGRangesFromDataFrame(hg38_rg_editnoY,
                                          ignore.strand = TRUE,
                                          keep.extra.columns = TRUE)

  data("copykit_obj_rle")
  copykit_data_proc <- as.data.frame(do.call(cbind,
                                             lapply(copykit_obj_rle,
                                                    inverse.rle)))
  copykit_obj <- CopyKit(assays = list(segment_ratios = copykit_data_proc),
                         rowRanges = rg_hg38)
  copykit_obj <- logNorm(copykit_obj)
  metadata(copykit_obj)$genome <- "hg38"
  colData(copykit_obj)$sample <- names(copykit_data_proc)
  return(copykit_obj)
}
