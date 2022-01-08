#' countBreakpoints
#'
#' Considers changes in the segment ratios as breakpoints.
#' Counts the breakpoints for each chromosome arm separately.
#'
#' @param scCNA
#'
#' @return The scCNA object with a column of breakpoint counts added to colData.
#' @export
#'
#' @keywords internal
#'
#' @importFrom dplyr pull bind_rows mutate select
#' @importFrom SummarizedExperiment rowRanges seqnames
#'
#' @examples
#' copykit_obj <- copykit_example_filtered()
#' copykit_obj <- .countBreakpoints(copykit_obj)
.countBreakpoints <- function(scCNA) {

    # bindings for NSE
    arm <- chrarm <- NULL

    rg_chr <- SummarizedExperiment::rowRanges(scCNA) %>%
        as.data.frame() %>%
        dplyr::mutate(chrarm = paste0(seqnames, arm)) %>%
        dplyr::select(chrarm)

    dat_seg_cp <- segment_ratios(scCNA)

    # split by chrom
    message("Counting breakpoints.")
    dat_seg_split <- split(dat_seg_cp, dplyr::pull(rg_chr, chrarm))

    brkpt_by_chrom <-
        lapply(dat_seg_split, function(x) {
            apply(x, 2, function(i) {
                length(rle(i)$values) - 1
            }) %>%
                unlist()
        })

    brkpt_by_chrom_df <- dplyr::bind_rows(brkpt_by_chrom) %>%
        t() %>%
        as.data.frame()

    brkpt_count <- rowSums(brkpt_by_chrom_df)

    # making sure order is identical
    brkpt_count <- brkpt_count[SummarizedExperiment::colData(scCNA)$sample]

    SummarizedExperiment::colData(scCNA)$breakpoint_count <-
        brkpt_count

    return(scCNA)
}
