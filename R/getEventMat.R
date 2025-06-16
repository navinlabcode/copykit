#' getEventMat
#'
#' Merge aconsecutive genomic bins that have identical copy number profiles across all subclonal populations
#' If the inter-segmental distance between two consecutively merged segments falls within the threshold (default=2), 
#' the middle bin is treated as the ‘common breakpoint’.
#'
#' @param scCNA A scCNA object.
#' @param bin_adj An integer with the threshold below which 2 breakpoints will be merged
#' @param ploidy_trunc An integer with the threshold above which the copy number will be set to this number
#' @details \code{getEventMat} uses consensus CN matrix to create a event matrix
#' Merge aconsecutive genomic bins that have identical copy number profiles across all subclonal populations
#' If the inter-segmental distance between two consecutively merged segments falls within the threshold (default=2),
#' the middle bin is treated as the ‘common breakpoint’.
#'
#' @return A copy number event matrix
#' @export
#'
#'
#' @importFrom dplyr as_tibble rename mutate select
#' @importFrom tibble remove_rownames
#' @importFrom SummarizedExperiment rowRanges seqnames
#'
#' @examples
#' set.seed(1000)
#' copykit_obj <- copykit_example_filtered()[, sample(40)]
#' event_mat <- getEventMat(copykit_obj)
getEventMat <- function(
                        scCNA,          # consensus CN matrix of which will be converted to event matrix
                        bin_adj = 2,    # number of bins allowed to be adjusted to consider as the same breakpoint
                        ploidy_trunc = 8   # maximum integer value, all integer value larger than this will be set to this
                        ){
  
  ## trunc integer matrix
  seg_df = copykit::consensus(scCNA)
  seg_df[seg_df>=ploidy_trunc] = ploidy_trunc
  
  
  intmat <- SummarizedExperiment::rowRanges(scCNA) %>%
    dplyr::as_tibble() %>%
    dplyr::select(seqnames, start, end, arm) %>%
    cbind(seg_df) %>%
    tibble::remove_rownames()
  
  ## merge segments
  res_int <- as.data.frame(intmat[1,])
  for(i in 2:nrow(intmat)){
    if(identical(as.character(intmat[i,-c(1:4)]), as.character(intmat[i-1,-c(1:4)]))){
      next
    }else{
      res_int <- rbind(res_int, intmat[i,])
    }
  }
  res_int$bin <- as.numeric(rownames(res_int))
  
  ## finding common breakpoints
  res_int_cbp <- as.data.frame(res_int[1,])
  for(i in 2:(nrow(res_int)-1)){
    if(res_int$bin[i+1]-res_int$bin[i]<=bin_adj){
      next
    }else{
      res_int_cbp <- rbind(res_int_cbp, res_int[i,])
    }
  }
  
  res_int_cbp$n.bins=c(res_int_cbp$bin[-1], nrow(scCNA)+1) - res_int_cbp$bin
  
  res_int_cbp$end.pos = intmat$end[as.numeric(res_int_cbp$bin)+res_int_cbp$n.bins-1]
  res_int_cbp$end.chr = intmat$seqnames[as.numeric(res_int_cbp$bin)+res_int_cbp$n.bins-1]
  res_df <- res_int_cbp %>%
    dplyr::mutate(chrom.arm = if_else( (intmat$arm[bin+n.bins-1]==arm) & (intmat$seqnames[bin+n.bins-1]==seqnames), arm, "-")) %>%
    dplyr::rename(start.chr=seqnames, 
                  start.pos=start) %>%
    dplyr::select(starts_with("start"), end.chr, end.pos, bin, n.bins, chrom.arm, everything(), -arm, -end) %>%
    tibble::remove_rownames()
  
  return(res_df)
}
