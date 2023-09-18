#' calcInteger()
#'
#' Calculates the integer copy number profile for each single cell
#'
#' @param scCNA The CopyKit object.
#' @param assay String with the name of the assay to pull data from to calculate
#' integers.
#' @param method Method used to scale the ratio values to integer.
#' @param ploidy_value If method of choice is 'fixed' a ploidy value should be
#' provided.
#' @param name String specifying the name to be used to store the result in the
#' reducedDims of the output.
#' @param penalty An integer passed on to scquantum::ploidy.inference()
#' penalty argument
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying how the function
#' should be parallelized.
#'
#' @details
#'
#' CopyKit support the following methods for calculating integer copy number
#' matrices
#' \itemize{
#' \item{fixed:} When method argument is set to 'fixed' copykit extracts the
#' segment means from the scCNA object and multiplies those means by the value
#' provided in the argument ploidy_value.
#'
#'
#' \item{scquantum:} When the method argument is set to 'scquantum', CopyKit
#' applies \code{\link[scquantum]{ploidy.inference}} function to perform a
#' sample wise calculation returning the estimated compuational ploidy for
#' every single cell
#' }
#'
#' @return The CopyKit object with an assay slot named 'integer' that contains
#' a data frame with cells as columns and integerized segments as rows. And, in
#' case of method = 'scquantum' CopyKit adds three new elements to \code{colData}
#'  named 'ploidy' and 'ploidy_score' and the 'confidence ratio' obtained from
#'  scquantum for each cell.
#'
#' @export
#'
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment assay assayNames colData rowRanges
#' @importFrom scquantum ploidy.inference timeseries.iod
#'
#' @examples
#' copykit_obj <- mock_bincounts(ncells_diploid = 0, ncells = 10)
#' copykit_obj <- calcInteger(copykit_obj, method = "scquantum")

calcInteger <- function(scCNA,
                        assay = c("bincounts",
                                  'smoothed_bincounts',
                                  "segment_ratios"),
                        method = "fixed",
                        ploidy_value = NULL,
                        name = "integer",
                        penalty = 25,
                        BPPARAM = bpparam()) {
  # args
  assay = match.arg(assay)

  if ('smoothed_bincounts' %in% assayNames(scCNA) && assay == 'bincounts'
      && method == 'scquantum') {
    warning("CopyKit detected that knnSmooth() has been performed.")
    warning("If working with knnSmooth datasets we recommend using the assay 'smoothed_bincounts'")
  }

  # getting datasets
  if (assay == 'bincounts') {
    bin <- SummarizedExperiment::assay(scCNA, 'bincounts')
  }

  if (assay %in% c('smoothed_bincounts','segment_ratios')) {
    bin <- SummarizedExperiment::assay(scCNA, 'smoothed_bincounts')
  }

  seg <- SummarizedExperiment::assay(scCNA, 'segment_ratios')

  if (!is.null(ploidy_value)) {
    if (method == "fixed") {
      if (is.null(ploidy_value) && !is.numeric(ploidy_value)) {
        stop("Method fixed requires a numeric value for ploidy_value.")
      }

      message("Scaling ratio values by ploidy value ",
              ploidy_value)

      # ploidy values are added to colData information
      SummarizedExperiment::colData(scCNA)$ploidy <-
        ploidy_value

      # saving ploidy scaling method
      S4Vectors::metadata(scCNA)$ploidy_method <- "fixed"
    }
  }

  if (method == "metadata") {
    # method metadata just allows the segment ratios to be integerized based
    # on the values for each cell in the colData(scCNA)$ploidy information.
    if (!is.null(colData(scCNA)$ploidy)) {
      message('Calculating integer values based on colData(scCNA)$ploidy info.')
    } else {
      stop("Method 'metadata' requires colData(scCNA)$ploidy information.")
    }

  }

  # logic for scquantum method
  if (method == "scquantum") {
    rg <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))

    # method scquantum with input bincounts
    if (assay %in% c('bincounts', 'smoothed_bincounts')) {
      sc_quants <-
        BiocParallel::bplapply(
          assay(scCNA, assay),
          scquantum::ploidy.inference,
          chrom = rg$seqnames,
          start = rg$start,
          end = rg$end,
          penalty = penalty,
          BPPARAM = BPPARAM
        )
    }

    # method scquantum with input segment ratios
    if (assay == 'segment_ratios') {

      sc_quants <- BiocParallel::bplapply(seq_along(seg), function(z) {
        # extracting segments rle id and lengths
        segnums <- cumsum(c(TRUE, abs(diff(seg[, z])) > 0.00001))
        seg_length <- rle(seg[, z])$lengths

        # extracting segment-wise means and index of dispersion
        seg_bins_mean <- tapply(bin[, z], segnums, mean)

        if (any(seg_length <= 3)) {
          iod.est <- scquantum::timeseries.iod(bin[,z])
        } else {
          iod.est <-
            tapply(bin[, z], segnums, scquantum::timeseries.iod)
        }

        # bincount mean estimate
        mean.est <- mean(bin[, z])

        estimates <- scquantum::ploidy.inference(
          x = seg_bins_mean,
          chrom = NULL,
          start = NULL,
          end = NULL,
          seg_length = seg_length,
          iod = iod.est,
          mean_bincount = mean.est,
          do_segmentation = FALSE
        )

      })


    }

    # extracting ploidies from sc_quantum object
    sc_ploidies <-
      vapply(sc_quants, function(x)
        x$ploidy, numeric(1))

    # extracting ploidies from sc_quantum object
    sc_confidence <-
      vapply(sc_quants, function(x)
        x$confidence_ratio, numeric(1))

    # calculating ploidy score from scquantum confidence ratio
    ploidy_score <- abs(1-sc_confidence)

    SummarizedExperiment::colData(scCNA)$ploidy <- sc_ploidies
    SummarizedExperiment::colData(scCNA)$confidence_ratio <-
      sc_confidence
    SummarizedExperiment::colData(scCNA)$ploidy_score <-
      ploidy_score

  }

  # check to guarantee multiplication
  if (!identical(names(bin), colData(scCNA)$sample)) {
    stop("Order of cells in segment_ratios and colData() is not identical.")
  }

  # obtain the matrix of integer values by multiplying the seg ratios
  # by the diagonal of the ploidy colData vector
  int_values <-
    round(as.matrix(seg) %*% diag(colData(scCNA)$ploidy)) %>%
    as.data.frame()

  # recovering names
  names(int_values) <- names(seg)

  SummarizedExperiment::assay(scCNA, name) <- int_values

  return(scCNA)
}
