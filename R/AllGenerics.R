###################################################################
# Getters/setters for scCNA

#' @export
setGeneric("segment_ratios", function(x, ...)
  standardGeneric("segment_ratios"))

#' @export
#' @rdname CopyKit-class
setGeneric("ratios", function(x, ...)
  standardGeneric("ratios"))

#' @export
#' @rdname CopyKit-class
setGeneric("bincounts", function(x, ...)
  standardGeneric("bincounts"))

#' @export
setGeneric("consensus", function(x, ...)
  standardGeneric("consensus"))

#' @export
setGeneric("consensus<-", function(x, ..., value)
  standardGeneric("consensus<-"))

#' @export
setGeneric("phylo", function(x, ...)
  standardGeneric("phylo"))

#' @export
setGeneric("phylo<-", function(x, ..., value)
  standardGeneric("phylo<-"))

#' @export
setGeneric("consensusPhylo", function(x, ...)
  standardGeneric("consensusPhylo"))

#' @export
setGeneric("consensusPhylo<-", function(x, ..., value)
  standardGeneric("consensusPhylo<-"))

#' @export
setGeneric("distMat", function(x, ...)
  standardGeneric("distMat"))

#' @export
setGeneric("distMat<-", function(x, ..., value)
  standardGeneric("distMat<-"))

#' @export
setGeneric("graph", function(x, ...)
  standardGeneric("graph"))

#' @export
setGeneric("graph<-", function(x, ..., value)
  standardGeneric("graph<-"))
