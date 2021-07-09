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
setMethod("consensusPhylo", "scCNA", function(x) {
  # accessor for the consensusPhylo slot
  out <- x@consensusPhylo
  out
})

#' @export
setReplaceMethod("consensusPhylo", "scCNA", function(x, value) {
  # setter method for consensusPhylo slot
  x@consensusPhylo <- value
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
#' @importFrom ape Ntip Nnode
#' @importMethodsFrom SingleCellExperiment show
setMethod("show", "scCNA", function(object) {
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

#' @export
# combine fun needed for Merge levels
.combine.func <-
  function(diff,
           vecObs,
           vecPredNow,
           mnNow,
           mn1,
           mn2,
           pv.thres = 0.0001,
           thresAbs = 0)
  {
    #observed values in the first segment
    vec1 <- vecObs[which(vecPredNow == mn1)]
    #observed values in the second segment
    vec2 <- vecObs[which(vecPredNow == mn2)]

    #if difference between segment medians does not exceed thresAbs, then set pv=1
    if (diff <= thresAbs) {
      pv <- 1
    }
    #otherwise test for difference in mean based on observed values
    else {
      if ((length(vec1) > 10 &
           length(vec2) > 10) |
          sum(length(vec1), length(vec2)) > 100) {
        pv <- wilcox.test(vec1, vec2)$p.value
      }
      else{
        pv <- wilcox.test(vec1, vec2, exact = T)$p.value
      }	#/10^max(mn1,mn2)
      if (length(vec1) <= 3 | length(vec2) <= 3) {
        pv <- 0
      }
    }
    index.merged <- numeric()
    #if p-value exceeds pv.thres
    if (pv > pv.thres) 	{
      #combine observed values
      vec <- c(vec1, vec2)
      # Index values to be updated
      index.merged <-
        which((vecPredNow == mn1) | (vecPredNow == mn2))
      #update predicted values by median of the observed values
      vecPredNow[index.merged] <- median(vec, na.rm = TRUE)
      #update segment medians  median of the observed values and remove one of the duplicates
      mnNow[which((mnNow == mn1) |
                    (mnNow == mn2))] <- median(vec, na.rm = TRUE)
      mnNow <- unique(mnNow)
    }
    list(mnNow = mnNow,
         vecPredNow = vecPredNow,
         pv = pv)
  }

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

#' @export
#' @keywords internal
# Merge Levels function
.MergeLevels <-
  function(vecObs,
           vecPred,
           pv.thres = 0.0001,
           ansari.sign = 0.05) {
    # Initializing threshold vector for keeping track of thresholds
    sq <- numeric()

    #initializing threshold index (threshold count)
    j <- 0

    #initializing ansari p-values to keep track of ansari p-values for each threshold in sq
    ansari <- numeric()

    # Initialize levels count
    lv <- numeric()

    # Start with threshold 0.05, and flag=0 indicating significance not yet reached, backtracking not begun
    thresAbs <- 0.05
    flag <- 0
    while (1) {
      j <- j + 1
      # Save current threshold
      sq[j] <- thresAbs

      # temporary predicted values (to be updated)
      vecPredNow <- vecPred

      #unmissing unique segment medians
      mnNow <- unique(vecPred)
      mnNow <- mnNow[!is.na(mnNow)]

      #continuing indicator otherwise get out of the loop
      cont <- 0

      while (cont == 0 & length(mnNow) > 1) {
        mnNow <- sort(mnNow)  #currennt sorted vector of means
        n <- length(mnNow)  # number of means in mnNow
        # cat("\r",n,":",length(unique(vecPred)))
        # Get distances translated to copy number differences
        # Only distances to closest levels
        d <- (2 * 2 ^ mnNow)[-n] - (2 * 2 ^ mnNow)[-1]

        #order distance between means with the closest on top and corresponding indices
        dst <-
          cbind(abs(d)[order(abs(d))], (2:n)[order(abs(d))], (1:(n - 1))[order(abs(d))])

        #for each pair of means
        for (i in 1:nrow(dst)) 	{
          #set continuity index to "NOT continue" (=1)
          cont <- 1
          #test for combining of the two segment means
          out <-
            .combine.func(
              diff = dst[i, 1],
              vecObs,
              vecPredNow,
              mnNow,
              mn1 = mnNow[dst[i, 2]],
              mn2 = mnNow[dst[i, 3]],
              pv.thres = pv.thres,
              thresAbs = thresAbs
            )
          #if combine?
          if (out$pv > pv.thres) {
            #set continuity index to "YES" (=0) and break out of the current pairs loop
            cont <- 0

            #update predicted values and segments
            vecPredNow <- out$vecPredNow
            mnNow <- out$mnNow
            break
          }
        }
      }
      # When done merging for a given threshold, test for significance
      ansari[j] <-
        ansari.test(sort(vecObs - vecPredNow), sort(vecObs - vecPred))$p.value
      if (is.na(ansari[j])) {
        ansari[j] <- 0
      }
      lv[j] <- length(mnNow) # get number of levels

      # If p.value is less than the significance threshold, set backtracking flag=1 (backtracking on)
      if (ansari[j] < ansari.sign) {
        flag <- 1
      }
      if (flag == 2) {
        break
      }

      # If backtracking is on, a smaller threshold is attempted
      if (flag) {
        # If backtracking is on and p.value is higher than sign threshold, stop
        if (ansari[j] > ansari.sign | thresAbs == 0) {
          # Don't merge at all if all tested threshold including 0 is significant
          if (ansari[j] <= ansari.sign) {
            vecPredNow <- vecPred
            mnNow <- unique(vecPred)
            mnNow <- mnNow[!is.na(mnNow)]
          }
          break
        }

        # Attempt smaller threshold
        else {
          thresAbs <- signif(thresAbs - 0.005, 3)
        }
      }
      else {
        thresAbs <-
          thresAbs + 0.1
      } # Increase threshold if backtracking is not on

      # Control step so function won't keep running, max threshold = 1 and if sign not reached, threshold = 0
      if (thresAbs >= 1) {
        thresAbs <- 0
        flag <- 2
      }
    }


    # Return list of results
    return(list(
      vecMerged = vecPredNow,
      mnNow = mnNow,
      sq = sq,
      ansari = ansari
    ))
  }

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
    c(
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
    names = paste0('c', 1:51)
  )
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
parCor <- function(x, BPPARAM=BiocParallel::bpparam()){
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
  
  return(res_parcor_reserve)
}
