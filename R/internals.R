#' @title CopyKit internal functions.
#' @name CopyKit-internals
#'
#' @description This document establish setters and getters to facilitate access
#' to fields for the CopyKit class object. The functions provided here are
#' in addition to setters and getters available from the SingleCellExperiment
#' class
#'
#' @section Getters:
#' \describe{
#' \item{\code{segment_ratios}:}{Returns a data frame of normalized segment
#' ratio means.}
#' \item{\code{ratios}:}{Returns a data frame of normalized ratio means.}
#' \item{\code{bincounts}:}{Returns a data frame of binned bincounts.}
#' \item{\code{consensus}:}{Returns a data frame of normalized segment
#' ratio means for the consensus matrix.}
#' \item{\code{phylo}:}{Returns a phylo class object with a phylogenetic tree.}
#' \item{\code{consensusPhylo}:}{Returns a phylo class object with a
#' phylogenetic tree from the consensus matrix.}
#' }
#'
#' @rdname internals
NULL

# -------------------
# global variables
# -------------------
utils::globalVariables(c('hg38_grangeslist', 'copykit_obj_filt_rle',
                         'copykit_obj_filt_umap', 'hg19_genes', 'hg38_genes',
                         'hg19_rg'))

#' @export
#' @rdname internals
#' @name segment_ratios
#' @aliases segment_ratios,CopyKit-method
#' @param x CopyKit object.
setMethod("segment_ratios", "CopyKit", function(x, withDimnames = TRUE) {
    # accessor for the segment_ratios data within the assay slot
    SummarizedExperiment::assay(x, "segment_ratios")
})


#' @export
#' @rdname internals
#' @name ratios
#' @aliases ratios,CopyKit-method
#' @param x CopyKit object.
setMethod("ratios", "CopyKit", function(x, withDimnames = TRUE) {
    # accessor for the ratios data slot
    SummarizedExperiment::assay(x, "ratios")
})

#' @export
#' @rdname internals
#' @name bincounts
#' @aliases bincounts,CopyKit-method
#' @param x CopyKit object.
setMethod("bincounts", "CopyKit", function(x, withDimnames = TRUE) {
    # accessor for the bincounts data slot
    SummarizedExperiment::assay(x, "bincounts")
})

#' @export
#' @rdname internals
#' @name consensus
#' @aliases consensus,CopyKit-method
#' @param x CopyKit object.
setMethod("consensus", "CopyKit", function(x, withDimnames = TRUE) {
    # accessor for the consensus data slot
    out <- x@consensus
    out
})


#' @export
#' @rdname internals
#' @name consensus<-
#' @aliases consensus<-,CopyKit-method
#' @param x CopyKit object.
setReplaceMethod("consensus", "CopyKit", function(x, value) {
    # setter method for phylo slot
    x@consensus <- value
    x
})

#' @export
#' @rdname internals
#' @name phylo
#' @aliases phylo,CopyKit-method
#' @param x CopyKit object.
setMethod("phylo", "CopyKit", function(x) {
    # accessor for the phylo slot
    out <- x@phylo
    out
})

#' @export
#' @rdname internals
#' @name phylo<-
#' @aliases phylo<-,CopyKit-method
#' @param x CopyKit object.
setReplaceMethod("phylo", "CopyKit", function(x, value) {
    # setter method for phylo slot
    x@phylo <- value
    x
})

#' @export
#' @rdname internals
#' @name consensusPhylo
#' @aliases consensusPhylo,CopyKit-method
#' @param x CopyKit object.
setMethod("consensusPhylo", "CopyKit", function(x) {
    # accessor for the consensusPhylo slot
    out <- x@consensusPhylo
    out
})

#' @export
#' @rdname internals
#' @name consensusPhylo<-
#' @aliases consensusPhylo<-,CopyKit-method
#' @param x CopyKit object.
setReplaceMethod("consensusPhylo", "CopyKit", function(x, value) {
    # setter method for consensusPhylo slot
    x@consensusPhylo <- value
    x
})

#' @export
#' @rdname internals
#' @name distMat
#' @aliases distMat,CopyKit-method
#' @param x CopyKit object.
setMethod("distMat", "CopyKit", function(x) {
    # accessor for the distMat slot
    out <- x@distMat
    out
})

#' @export
#' @rdname internals
#' @name distMat<-
#' @aliases distMat<-,CopyKit-method
#' @param x CopyKit object.
setReplaceMethod("distMat", "CopyKit", function(x, value) {
    # setter method for distMat slot
    x@distMat <- value
    x
})

#' @export
#' @rdname internals
#' @name graph
#' @aliases graph,CopyKit-method
#' @param x CopyKit object.
setMethod("graph", "CopyKit", function(x) {
    # accessor for the getGraph slot
    out <- x@graph
    out
})

#' @export
#' @rdname internals
#' @name graph<-
#' @aliases graph<-,CopyKit-method
#' @param x CopyKit object.
setReplaceMethod("graph", "CopyKit", function(x, value) {
    # setter method for getGraph slot
    x@graph <- value
    x
})

#' @export
#' @rdname internals
#' @name show
#' @aliases show,CopyKit-method
#' @param x CopyKit object.
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

#' Not in
#'
#' @export
#' @rdname internals
#' @keywords internal
# thanks to https://peter.solymos.org/code/2016/11/26/
#how-to-write-and-document-special-functions-in-r.html
`%!in%` <- function(x, table) !(match(x, table, nomatch = 0) > 0)

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
.invft <- function(y) (y^2 - 1)^2 / (4 * y^2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# color palettes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# heatmaps ocean.balance
# The default color scale for Heatmap is from the package pals by Kevin Wright
# https://github.com/kwstat/pals (License GPL-3)
# The ocean.balance color palette (Licence MIT) was developed by Kristen Thyng
# http://dx.doi.org/10.5670/oceanog.2016.66
# Both awesome palette packages and you should check them out.

#' @name ocean_balance_hex
#' @aliases ocean_balance_hex
#' @export
#' @docType methods
#' @return Hexadecimal values for ocean.balance function
#' @rdname internals
ocean_balance_hex <- c(
    "#181C43",
    "#21295F",
    "#27347D",
    "#28409F",
    "#1950BA",
    "#0A65BD",
    "#1F77BB",
    "#3787BA",
    "#4F96BA",
    "#6AA4BC",
    "#87B2C1",
    "#A3BFC9",
    "#BECDD3",
    "#D8DCDE",
    "#F0ECEB",
    "#E8D7D2",
    "#E1C1B8",
    "#DAAD9E",
    "#D49884",
    "#CD8369",
    "#C66E51",
    "#BF583B",
    "#B7412A",
    "#AB2924",
    "#991527",
    "#830E29",
    "#6B0E25",
    "#530D1D",
    "#3C0912"
)

#' @name ocean.balance
#' @aliases ocean.balance
#' @export
#' @docType methods
#' @return The default colors for heatmaps
#' @importFrom grDevices colorRampPalette
#' @rdname internals
ocean.balance <- function(n) {
    colorRampPalette(ocean_balance_hex)(n)
}

#' @name superclones_pal
#' @aliases superclones_pal
#' @export
#' @docType methods
#' @return a named vector of default colors for CopyKit superclones.
#' @rdname internals
# superclones palette
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
        names = paste0("s", 1:24)
    )
}

#' @name subclones_pal
#' @aliases subclones_pal
#' @export
#' @docType methods
#' @return a named vector of default colors for CopyKit subclones.
#' @rdname internals
# subclones palette
subclones_pal <- function() {
    structure(
        c(
            "gray",
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
        names = paste0("c", 0:51)
    )
}

#' @name find_scaffold_genes
#' @aliases find_scaffold_genes
#' @export
#' @docType methods
#' @return A data frame with the gene HUGO gene symbol and the position on the
#' relevant scaffold from the varbin pipeline.
#' @rdname internals
#' @keywords internal
find_scaffold_genes <- function(scCNA,
                                genes) {

    # lazydata bindings and NSE
    symbol <- gene <- NULL

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
    genes_features <- BiocGenerics::subset(
        genes_assembly,
        symbol %in% genes
    )

    all_genes <- genes_assembly$symbol %>%
        unlist() %>%
        unname()

    missing_genes <- genes[genes %!in% all_genes]

    if (length(missing_genes) != 0) {
        warning(
            base::paste(
                "Genes:",
                paste(missing_genes,
                    collapse = ", "
                ),
                ",could not be found. Maybe you need to use a different gene alias?"
            )
        )
    }

    # finding overlaps
    olaps <-
        suppressWarnings(GenomicRanges::findOverlaps(genes_features,
            ranges,
            ignore.strand = TRUE
        ))

    # creating a data frame that will contain the genes and positions
    # (index) in the pipeline ranges.
    # some genes might overlap more than one range (more than one bin),
    # in this case only one will be kept
    df <-
        data.frame(
            gene = as.character(genes_features$symbol[S4Vectors::queryHits(olaps)]),
            pos = S4Vectors::subjectHits(olaps)
        ) %>%
        dplyr::distinct(gene, .keep_all = TRUE)

    # checking for genes that might have been excluded from the varbin pipeline
    blk_list <- genes[genes %!in% missing_genes]
    blk_list <- blk_list[blk_list %!in% df$gene]

    if (length(blk_list) != 0) {
        warning(
            base::paste(
                "Genes:",
                paste(blk_list,
                    collapse = ", "
                ),
                "are in excluded regions of the Varbin pipeline and can't be plotted."
            )
        )
    }

    return(df)
}


#' @author Alexander Davis
#' @name l2e.normal.sd
#' @aliases l2e.normal.sd
#' @export
#' @docType methods
#' @return A numeric vector with least squares sd estimation.
#' @rdname internals
l2e.normal.sd <- function(xs) {
    # Need at least two values to get a standard deviation
    stopifnot(length(xs) >= 2)
    optim.result <- stats::optimize(
        # L2E loss function
        f = function(sd) {
              # "Data part", the sample average of the likelihood
              -2 * mean(stats::dnorm(xs, sd = sd)) +
                  # "Theta part", the integral of the squared density
                  1 / (2 * sqrt(pi) * sd)
          },
        # Parameter: standard deviation of the normal distribution fit
        interval = c(0, diff(range(xs)))
    )
    return(optim.result$minimum)
}

#' @author Alexander Davis
#' @name overdispersion
#' @aliases overdispersion
#' @export
#' @docType methods
#' @return A numerci vector with the estimation of the index of dispersion,
#' which is used when estimating standard errors for each segment mean
#' @rdname internals
overdispersion <- function(v) {
    # 3 elements, 2 differences, can find a standard deviation
    stopifnot(length(v) >= 3)
    # Differences between pairs of values
    y <- v[-1]
    x <- v[-length(v)]
    # Normalize the differences using the sum. The result should be around zero,
    # plus or minus square root of the index of dispersion
    vals.unfiltered <- (y - x) / sqrt(y + x)
    # Remove divide by zero cases, and--considering this is supposed to be count
    # data--divide by almost-zero cases
    vals <- vals.unfiltered[y + x >= 1]
    # Check that there's anything left
    stopifnot(length(vals) >= 2)
    # Assuming most of the normalized differences follow a normal distribution,
    # estimate the standard deviation
    val.sd <- l2e.normal.sd(vals)
    # Square this standard deviation to obtain an estimate of the index of
    # dispersion
    iod <- val.sd^2
    # subtract one to get the overdispersion criteria
    iod.over <- iod - 1
    # normalizing by mean bincounts
    iod.norm <- iod.over / mean(v)
    return(iod.norm)
}

#' @author Junke Wang
#' @name parCor
#' @aliases parCor
#' @export
#' @docType methods
#' @return A matrix with the pairwise correlation from the segment ratio means.
#' @rdname internals
parCor <- function(x, BPPARAM = BiocParallel::bpparam()) {
    ncol <- ncol(x)

    ## skip parallelization if # of cell less than 2000
    if (ncol < 2001) {
        return(cor(x))
    }

    nSplit <- floor(ncol / 1000)
    lSplit <- floor(ncol / nSplit)
    iSplit <- vector("list", nSplit)
    for (i in seq_len((nSplit - 1))) {
        iSplit[[i]] <- (lSplit * (i - 1) + 1):(lSplit * i)
    }
    iSplit[[nSplit]] <- (lSplit * (nSplit - 1) + 1):ncol

    comb <- expand.grid(1:nSplit, 1:nSplit)
    comb <- unique(t(apply(comb, 1, sort)))

    ## result matrix
    result <- BiocParallel::bplapply(
        X = 1:nrow(comb),
        FUN = function(i) {
            a <- list()
            a[[1]] <- iSplit[[comb[i, 1]]]
            a[[2]] <- iSplit[[comb[i, 2]]]
            a[[3]] <- cor(x[, a[[1]]], x[, a[[2]]])
            a
        },
        BPPARAM = BPPARAM
    )

    res_parcor_reserve <- matrix(, nrow = ncol, ncol = ncol)
    for (i in 1:length(result)) {
        res_parcor_reserve[result[[i]][[1]], result[[i]][[2]]] <- result[[i]][[3]]
        res_parcor_reserve[result[[i]][[2]], result[[i]][[1]]] <- t(result[[i]][[3]])
    }

    colnames(res_parcor_reserve) <- colnames(x)
    rownames(res_parcor_reserve) <- colnames(x)
    return(res_parcor_reserve)
}

#' CopyKit Example
#'
#' @name copykit_example
#' @aliases copykit_example
#' @docType data
#' @rdname data
#' @return A CopyKit object with data from BL1 sample from the CopyKit
#' manuscript.
#' @details This function is largely used for examples of CopyKit functions.
#'
#' This function returns a CopyKit object with the unfiltered dataset
#' of the sample BL1. This sample was processed with ACT and is a liver
#' metastasis from a primary breast cancer patient.
#'
#' The CopyKit object contain only the segment ratio mean of each single cell
#' and, therefore, is not compatible with the functions from the runVarbin
#' module.
#' @export
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment colData<-
#' @keywords internal
copykit_example <- function() {

    # bindings for NSE
    chr <- NULL

    # ranges
    hg38_rg <- hg38_grangeslist[["hg38_200kb"]]
    hg38_rg <- subset(hg38_rg, seqnames != "chrY")

    copykit_obj_rle <- copykit_obj_rle
    copykit_data_proc <- as.data.frame(do.call(
        cbind,
        lapply(
            copykit_obj_rle,
            inverse.rle
        )
    ))
    copykit_obj <- CopyKit(
        assays = list(segment_ratios = copykit_data_proc),
        rowRanges = hg38_rg
    )
    copykit_obj <- logNorm(copykit_obj)
    metadata(copykit_obj)$genome <- "hg38"
    colData(copykit_obj)$sample <- names(copykit_data_proc)
    return(copykit_obj)
}

#' CopyKit Example
#'
#' @name copykit_example_filtered
#' @aliases copykit_example_filtered
#' @docType data
#' @rdname data
#' @return A CopyKit object with data from BL1 sample from the CopyKit
#' manuscript with diploid and noise cells removed.
#' @details This function is largely used for examples of CopyKit functions.
#'
#' This function returns a CopyKit object with the unfiltered dataset
#' of the sample BL1. This sample was processed with ACT and is a liver
#' metastasis from a primary breast cancer patient.
#'
#' The CopyKit object contain only the segment ratio mean of each single cell
#' and, therefore, is not compatible with the functions from the runVarbin
#' module.
#' @export
#' @importFrom S4Vectors metadata<-
#' @importFrom SummarizedExperiment colData<-
#' @keywords internal
#' @importFrom S4Vectors metadata
#' @importFrom SingleCellExperiment reducedDim<-
copykit_example_filtered <- function() {
    # NSE bindings
    chr <- NULL

    # ranges

    hg38_rg <- hg38_grangeslist[["hg38_200kb"]]
    hg38_rg <- subset(hg38_rg, seqnames != "chrY")

    # lazydata loaded objects
    copykit_obj_filtered <- copykit_obj_filt_rle
    copykit_obj_filtered_umap <- copykit_obj_filt_umap

    copykit_data_proc_filt <- as.data.frame(do.call(
        cbind,
        lapply(
            copykit_obj_filtered,
            inverse.rle
        )
    ))
    copykit_obj_filt <- CopyKit(
        assays = list(segment_ratios = copykit_data_proc_filt),
        rowRanges = hg38_rg
    )
    copykit_obj_filt <- logNorm(copykit_obj_filt)
    metadata(copykit_obj_filt)$genome <- "hg38"
    metadata(copykit_obj_filt)$suggestedK <- 10
    colData(copykit_obj_filt)$sample <- names(copykit_data_proc_filt)
    reducedDim(copykit_obj_filt, "umap") <- copykit_obj_filtered_umap
    return(copykit_obj_filt)
}

#' mock_bincounts
#'
#' @name mock_bincounts
#' @param ncells A numeric with the total number of cells to simulate.
#' @param ncells_diploid A numeric with the number of diploid cells to simulate
#' @param position_gain A vector with the index of the bin counts in which
#' chromosomal gains will be added.
#' @param position_del A vector with the index of the bin counts in which
#' chromosomal deletions will be placed.
#' @param genome The assembly genome information to add to the metadata.
#' @param resolution The resolution of the scaffold to add to the object
#' metadata
#' @aliases mock_bincounts
#' @export
#' @docType data
#' @return A CopyKit object with simulated bincounts
#' @rdname data
#' @keywords internal
#' @importFrom stats rpois runif
mock_bincounts <- function(ncells = 30,
                           ncells_diploid = 5,
                           position_gain = 4900:5493,
                           position_del = 6523:7056,
                           genome = "hg38",
                           resolution = "220kb",
                           run_vst = TRUE,
                           run_segmentation = TRUE,
                           run_lognorm = TRUE) {
    hg38_rg <- switch(resolution,
        "55kb" = hg38_grangeslist[["hg38_50kb"]],
        "110kb" = hg38_grangeslist[["hg38_100kb"]],
        "195kb" = hg38_grangeslist[["hg38_175kb"]],
        "220kb" = hg38_grangeslist[["hg38_200kb"]],
        "280kb" = hg38_grangeslist[["hg38_250kb"]],
        "500kb" = hg38_grangeslist[["hg38_500kb"]],
        "1Mb" = hg38_grangeslist[["hg38_1Mb"]],
        "2.8Mb" = hg38_grangeslist[["hg38_2Mb"]]
    )

    hg38_rg <- subset(hg38_rg, seqnames != "chrY")

    ncells <- ncells
    ncells_diploid <- ncells_diploid
    ncells_aneuploid <- ncells - ncells_diploid
    nbins <- length(hg38_rg)
    position_gain <- position_gain
    nbins_gain <- length(position_gain)
    position_del <- position_del
    nbins_del <- length(position_del)

    # creating mock diploid and aneuploid cell
    mock_diploid <- withr::with_seed(123, rpois(nbins, 50))
    mock_aneuploid <- withr::with_seed(123, rpois(nbins, 50))

    # adding events to aneuploid cells
    mock_aneuploid[position_gain] <- withr::with_seed(123,
                                                      rpois(nbins_gain, 100))
    mock_aneuploid[position_del] <- withr::with_seed(123, rpois(nbins_del, 25))

    # creating the cell counts matrix
    m <- matrix(c(mock_diploid, mock_aneuploid), ncol = 2)
    rep_pop <- c(rep(1, ncells_diploid), rep(2, ncells_aneuploid))
    mock_counts <- as.data.frame(m[, rep_pop])

    # adding some uniform error to avoid all cells having the same variance.
    mock_counts <- mock_counts + withr::with_seed(123,
                                                  runif(nbins * ncells, -5, 5))

    # creating copykit object with mock counts
    copykit_obj_bincounts <- CopyKit(
        assays = list(bincounts = mock_counts),
        rowRanges = hg38_rg
    )

    # adding metadata elements
    metadata(copykit_obj_bincounts)$genome <- genome
    metadata(copykit_obj_bincounts)$resolution <- resolution
    colData(copykit_obj_bincounts)$sample <- names(bincounts(copykit_obj_bincounts))
    colData(copykit_obj_bincounts)$ground_truth <- rep(
        c(FALSE, TRUE),
        c(
            ncells_diploid,
            ncells_aneuploid
        )
    )

    if (run_vst == TRUE) {
        copykit_obj_bincounts <- runVst(copykit_obj_bincounts)
    }

    if (run_segmentation == TRUE) {
        copykit_obj_bincounts <- runSegmentation(copykit_obj_bincounts)
    }

    if (run_lognorm == TRUE) {
        copykit_obj_bincounts <- logNorm(copykit_obj_bincounts)
    }

    return(copykit_obj_bincounts)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge Levels
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mergeLevels and combine.func are present in https://doi.org/10.1038/ng.3641
# They can also be found in the package aCGH from Peter Dimitrov
# https://www.bioconductor.org/packages/release/bioc/html/aCGH.html
# under GPL-2 licence

#' @name mergeLevels
#' @aliases mergeLevels
#' @export
#' @docType methods
#' @rdname internals
#' @keywords internal
#' @importFrom stats ansari.test wilcox.test
mergeLevels <- function (vecObs, vecPred, pv.thres = 1e-04, ansari.sign = 0.05,
                         thresMin = 0.05, thresMax = 0.5, verbose = 1, scale = TRUE)
{
    if (thresMin > thresMax) {
        cat("Error, thresMax should be equal to or larger than thresMin\n")
        return()
    }
    thresAbs = thresMin
    sq <- numeric()
    j = 0
    ansari = numeric()
    lv = numeric()
    flag = 0
    if (thresMin == thresMax) {
        flag = 2
    }
    else {
        l.step <- signif((thresMax - thresMin)/10, 1)
        s.step <- signif((thresMax - thresMin)/200, 1)
    }
    while (1) {
        if (verbose >= 1) {
            cat("\nCurrent thresAbs: ", thresAbs, "\n")
        }
        j = j + 1
        sq[j] <- thresAbs
        vecPredNow = vecPred
        mnNow = unique(vecPred)
        mnNow = mnNow[!is.na(mnNow)]
        cont = 0
        while (cont == 0 & length(mnNow) > 1) {
            mnNow = sort(mnNow)
            n <- length(mnNow)
            if (verbose >= 2) {
                cat("\r", n, ":", length(unique(vecPred)), "\t")
            }
            if (scale) {
                d <- (2 * 2^mnNow)[-n] - (2 * 2^mnNow)[-1]
            }
            else {
                d <- (mnNow)[-n] - (mnNow)[-1]
            }
            dst <- cbind(abs(d)[order(abs(d))], (2:n)[order(abs(d))],
                         (1:(n - 1))[order(abs(d))])
            for (i in 1:nrow(dst)) {
                cont = 1
                out = combine.func(diff = dst[i, 1], vecObs,
                                   vecPredNow, mnNow, mn1 = mnNow[dst[i, 2]],
                                   mn2 = mnNow[dst[i, 3]], pv.thres = pv.thres,
                                   thresAbs = if (scale) {
                                       2 * 2^thresAbs - 2
                                   }
                                   else {
                                       thresAbs
                                   })
                if (out$pv > pv.thres) {
                    cont = 0
                    vecPredNow = out$vecPredNow
                    mnNow = out$mnNow
                    break
                }
            }
        }
        ansari[j] = ansari.test(sort(vecObs - vecPredNow), sort(vecObs -
                                                                    vecPred))$p.value
        if (is.na(ansari[j])) {
            ansari[j] = 0
        }
        lv[j] = length(mnNow)
        if (flag == 2) {
            break
        }
        if (ansari[j] < ansari.sign) {
            flag = 1
        }
        if (flag) {
            if (ansari[j] > ansari.sign | thresAbs == thresMin) {
                break
            }
            else {
                thresAbs = signif(thresAbs - s.step, 3)
                if (thresAbs <= thresMin) {
                    thresAbs = thresMin
                }
            }
        }
        else {
            thresAbs = thresAbs + l.step
        }
        if (thresAbs >= thresMax) {
            thresAbs = thresMax
            flag = 2
        }
    }
    return(list(vecMerged = vecPredNow, mnNow = mnNow, sq = sq,
                ansari = ansari))
}

#' @name combine.func
#' @aliases combine.func
#' @export
#' @docType methods
#' @rdname internals
#' @keywords internal
combine.func <- function (diff, vecObs, vecPredNow, mnNow, mn1, mn2, pv.thres = 1e-04,
                          thresAbs = 0)
{
    vec1 = vecObs[which(vecPredNow == mn1)]
    vec2 = vecObs[which(vecPredNow == mn2)]
    if (diff <= thresAbs) {
        pv = 1
    }
    else {
        if ((length(vec1) > 10 & length(vec2) > 10) | sum(length(vec1),
                                                          length(vec2)) > 100) {
            pv = wilcox.test(vec1, vec2)$p.value
        }
        else {
            pv = wilcox.test(vec1, vec2, exact = TRUE)$p.value
        }
        if (length(vec1) <= 3 | length(vec2) <= 3) {
            pv = 0
        }
    }
    index.merged <- numeric()
    if (pv > pv.thres) {
        vec = c(vec1, vec2)
        index.merged = which((vecPredNow == mn1) | (vecPredNow ==
                                                        mn2))
        vecPredNow[index.merged] = median(vec, na.rm = TRUE)
        mnNow[which((mnNow == mn1) | (mnNow == mn2))] = median(vec,
                                                               na.rm = TRUE)
        mnNow = unique(mnNow)
    }
    list(mnNow = mnNow, vecPredNow = vecPredNow, pv = pv)
}

