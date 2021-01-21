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
setMethod("consensus", "scCNA", function(x, withDimnames = TRUE) {
  # accessor for the consensus data slot
  out <- x@consensus
  out
})

#' @export
setReplaceMethod("consensus", "scCNA", function(x, value) {
  # setter method for phylo slot
  x@consensus <- value
  x
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
    "Phylo:", phylo(object), "\n",
    "consensus: ", nrow(consensus(object)), " ", ncol(consensus(object))
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

#' @export
# combine fun needed for Merge levels
.combine.func=function(diff,vecObs, vecPredNow, mnNow, mn1, mn2, pv.thres=0.0001, thresAbs=0)
{
  #observed values in the first segment
  vec1=vecObs[which(vecPredNow==mn1)]
  #observed values in the second segment
  vec2=vecObs[which(vecPredNow==mn2)]

  #if difference between segment medians does not exceed thresAbs, then set pv=1
  if (diff<=thresAbs) {
    pv=1
  }
  #otherwise test for difference in mean based on observed values
  else {
    if((length(vec1) > 10 & length(vec2) > 10) | sum(length(vec1),length(vec2))>100){
      pv=wilcox.test(vec1,vec2)$p.value
    }
    else{pv=wilcox.test(vec1,vec2,exact=T)$p.value  }	#/10^max(mn1,mn2)
    if(length(vec1) <= 3 | length(vec2) <= 3){pv=0}
  }
  index.merged<-numeric()
  #if p-value exceeds pv.thres
  if (pv > pv.thres) 	{
    #combine observed values
    vec=c(vec1,vec2)
    # Index values to be updated
    index.merged=which((vecPredNow==mn1) | (vecPredNow==mn2))
    #update predicted values by median of the observed values
    vecPredNow[index.merged]=median(vec, na.rm=TRUE)
    #update segment medians  median of the observed values and remove one of the duplicates
    mnNow[which((mnNow==mn1) | (mnNow==mn2))]=median(vec, na.rm=TRUE)
    mnNow=unique(mnNow)
  }
  list(mnNow=mnNow, vecPredNow=vecPredNow, pv=pv)
}

#' @export
# Merge Levels function
.MergeLevels <- function(vecObs,vecPred,pv.thres=0.0001,ansari.sign=0.05){

  # Initializing threshold vector for keeping track of thresholds
  sq<-numeric()

  #initializing threshold index (threshold count)
  j=0

  #initializing ansari p-values to keep track of ansari p-values for each threshold in sq
  ansari=numeric()

  # Initialize levels count
  lv=numeric()

  # Start with threshold 0.05, and flag=0 indicating significance not yet reached, backtracking not begun
  thresAbs=0.05
  flag=0
  while (1){
    j=j+1
    # Save current threshold
    sq[j]<-thresAbs

    # temporary predicted values (to be updated)
    vecPredNow=vecPred

    #unmissing unique segment medians
    mnNow=unique(vecPred)
    mnNow=mnNow[!is.na(mnNow)]

    #continuing indicator otherwise get out of the loop
    cont=0

    while(cont==0 & length(mnNow)>1) {

      mnNow=sort(mnNow)  #currennt sorted vector of means
      n <- length(mnNow)  # number of means in mnNow
      # cat("\r",n,":",length(unique(vecPred)))
      # Get distances translated to copy number differences
      # Only distances to closest levels
      d<-(2*2^mnNow)[-n]-(2*2^mnNow)[-1]

      #order distance between means with the closest on top and corresponding indices
      dst<-cbind(abs(d)[order(abs(d))],(2:n)[order(abs(d))],(1:(n-1))[order(abs(d))])

      #for each pair of means
      for (i in 1:nrow(dst)) 	{
        #set continuity index to "NOT continue" (=1)
        cont=1
        #test for combining of the two segment means
        out=.combine.func(diff=dst[i,1],vecObs, vecPredNow, mnNow, mn1=mnNow[dst[i,2]], mn2=mnNow[dst[i,3]], pv.thres=pv.thres, thresAbs=thresAbs)
        #if combine?
        if (out$pv > pv.thres) {

          #set continuity index to "YES" (=0) and break out of the current pairs loop
          cont=0

          #update predicted values and segments
          vecPredNow=out$vecPredNow
          mnNow=out$mnNow
          break
        }
      }
    }
    # When done merging for a given threshold, test for significance
    ansari[j]=ansari.test(sort(vecObs-vecPredNow), sort(vecObs-vecPred))$p.value
    if(is.na(ansari[j])){ansari[j]=0}
    lv[j]=length(mnNow) # get number of levels

    # If p.value is less than the significance threshold, set backtracking flag=1 (backtracking on)
    if(ansari[j]<ansari.sign){
      flag=1
    }
    if(flag==2){ break }

    # If backtracking is on, a smaller threshold is attempted
    if (flag){

      # If backtracking is on and p.value is higher than sign threshold, stop
      if (ansari[j]>ansari.sign | thresAbs == 0){

        # Don't merge at all if all tested threshold including 0 is significant
        if (ansari[j] <= ansari.sign) {
          vecPredNow=vecPred
          mnNow=unique(vecPred)
          mnNow=mnNow[!is.na(mnNow)]
        }
        break
      }

      # Attempt smaller threshold
      else {thresAbs=signif(thresAbs-0.005,3) }
    }
    else {thresAbs=thresAbs+0.1} # Increase threshold if backtracking is not on

    # Control step so function won't keep running, max threshold = 1 and if sign not reached, threshold = 0
    if (thresAbs >= 1){
      thresAbs=0
      flag=2
    }
  }


  # Return list of results
  return(list(vecMerged=vecPredNow,mnNow=mnNow,sq=sq,ansari=ansari))
}
