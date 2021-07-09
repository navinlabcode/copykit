#' Plot phylogenetic tree
#'
#' Plots a phylogenetic tree.
#'
#' @author Darlan Conterno Minussi
#'
#' @param scCNA scCNA object.
#' @param label Character. Annotate tree tip by one of the element of metadata.
#' Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @param label_color List. Named list with colors for the label annotation.
#'  Must match label length
#' @param consensus Boolean. Indicates if the consensus phylogenetic tree should be plotted. Default to FALSE. The label/label_color option will be disabled.
#'
#' @return A phylogenetic tree visualization.
#'
#' @export
#'
#' @import ggtree
#' @examples
#'

plotPhylo <- function(scCNA,
					  label,
					  label_color,
					  consensus = FALSE){

  # consensus and not consensus logic
  if (consensus == FALSE) {
    tryCatch(
      phylo(scCNA),
      error = function(e) {
        message("No phylogeny detected in scCNA object.")
      },
      finally = {
        scCNA <- runPhylo(scCNA)
      }
    )

    tree <- phylo(scCNA)

  } else {

  	tryCatch(
      consensusPhylo(scCNA),
      error = function(e) {
        message("No consensus phylogeny detected in scCNA object.")
      },
      finally = {
        scCNA <- runConsensusPhylo(scCNA)
      }
    )

    tree <- consensusPhylo(scCNA)

  }


    # plotting
  if (is.null(label)) {
    # plotting without labels

    ggtree::ggtree(tree, laderize = F)

  } else {
  	# plotting with labels
  	if (length(label)>1) {
  		stop("Only one label can be visualized in the tree.")
  	}

    metadata <- SummarizedExperiment::colData(scCNA) %>%
      as.data.frame()

    if (!(label %in% colnames(metadata))) {
    	stop(paste0("Label ", label, " is not a column of the scCNA object."))
  	}

	metadata_anno_df <- metadata %>% dplyr::select(dplyr::all_of(label))

    if (consensus == FALSE) {
      
    	metadata_anno_df <- metadata_anno_df[tree$tip.label,]

    }


    if (consensus == TRUE) {
      # Uses the hidden consensus_by attribute from the calcConsensus function
      # to guarantee the same order
      cons_attr <- attr(consensus(scCNA), "consensus_by")

      metadata_anno_df <- metadata_anno_df[label] %>%
        dplyr::distinct()

      rownames(metadata_anno_df) <- metadata_anno_df %>%
        dplyr::pull(!!cons_attr)

      metadata_anno_df <-
        metadata_anno_df[names(consensus(scCNA)), , drop = FALSE]

    }

    if (is.null(label_color)) {
     
      #default colors superclones and subclones
      label_colors <- c(
        list(
          superclones = superclones_pal(),
          subclones = subclones_pal(),
          filtered = c("removed" = "#DA614D",
                       "kept" = "#5F917A"),
          is_normal = c("TRUE" = "#396DB3",
                        "FALSE" = "#11181D")
        ),
        label_colors
      )

        if (any(str_detect(
          label,
          c("superclones", "subclones", "filtered", "is_normal")
        ))) {
          # if label is one of the four above, uses the default specifed colors above
          label_color <- label_colors[[label]]

        } else if (is.numeric(dplyr::pull(metadata_anno_df, label)))  {
          # if label is a numeric vector
          n = 300
          min_v = min(dplyr::pull(metadata_anno_df, label))
          max_v = max(dplyr::pull(metadata_anno_df, label))

          label_color <-
            list(circlize::colorRamp2(
              seq(min_v, max_v, length = n),
              viridis::viridis(n, option = "D")
            ))
          names(label_color) <- label


        } else {
          # if label is not numeric

          elements <- metadata_anno_df %>%
            dplyr::pull(label) %>%
            unique() %>%
            as.character() %>%
            sort()

          n <- length(elements)
          hues <- seq(15, 375, length = n + 1)
          hex <- hcl(h = hues, l = 35, c = 100)[1:n]

          col <- structure(hex,
                           names = elements)

          label_color <- list(col)
          names(label_color) <- label


        }

      

    }


    list_samples <- split(rownames(metadata_anno_df), metadata_anno_df[,label])
    tree <- ggtree::groupOTU(tree, list_samples)
    p <- ggtree::ggtree(tree, ladderize=F) +
    	geom_tippoint(aes(color=group), size=1) +
    	scale_color_manual(values=label_color, name=label)

    if(consensus){
    	p<-p+geom_tiplab(aes(color=group), size=2)
    }

    p

}