#' Plot phylogenetic tree
#'
#' Plots a phylogenetic tree.
#'
#' @author Junke Wang
#'
#' @param scCNA scCNA object.
#' @param label Character. Annotate tree tip by one of the element of metadata.
#' Metadata can be accessed with \code{SummarizedExperiment::colData(scCNA)}
#' @param label_colors List. Named list with colors for the label annotation.
#'  Must match label length
#' @param consensus Boolean. Indicates if the consensus phylogenetic tree should be plotted. Default to FALSE. 
#'
#' @return A phylogenetic tree visualization.
#' @export
#'
#' @import ggtree
#' @examples
#'

plotPhylo <- function(scCNA,
					  label = NULL,
					  label_colors = NULL,
					  consensus = FALSE){

  # consensus and not consensus logic
  if (consensus == FALSE) {
    tryCatch(
      phylo(scCNA),
      error = function(e) {
        stop("No phylogeny detected in scCNA object. Please run runPhylo() first.")
      }
    )

    tree <- phylo(scCNA)

  } else {

  	tryCatch(
      consensusPhylo(scCNA),
      error = function(e) {
        stop("No consensus phylogeny detected in scCNA object. Please run runConsensusPhylo() first.")
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
  	if ( length(label)>1 ) {
  		stop("Only one label can be visualized in the tree.")
  	}

    metadata <- SummarizedExperiment::colData(scCNA) %>%
      as.data.frame()

    if (!(label %in% colnames(metadata))) {
    	stop(paste0("Label ", label, " is not a column of the scCNA object."))
  	}

	metadata_anno_df <- metadata %>% dplyr::select(dplyr::all_of(label))

    if (consensus == FALSE) {
      
    	metadata_anno_df <- metadata_anno_df[tree$tip.label,, drop=F]
    	size <- 1

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

      size <- 5

    }

    if (is.null(label_colors)) {
     
      #default colors superclones and subclones
      label_colors <- c(
        list(
          superclones = superclones_pal(),
          subclones = subclones_pal(),
          filtered = c("removed" = "#DA614D",
                       "kept" = "#5F917A"),
          is_normal = c("TRUE" = "#396DB3",
                        "FALSE" = "#11181D")
        )
      )

        if (any(str_detect(
          label,
          c("superclones", "subclones", "filtered", "is_normal")
        ))) {
          # if label is one of the four above, uses the default specifed colors above
          label_colors <- label_colors[[label]]

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

          label_colors <- col


        }

      

    } else if (is.list(label_colors)) {
    	label_colors <- label_colors[[label]]
    }


    list_samples <- split(rownames(metadata_anno_df), metadata_anno_df[,label])
    tree <- ggtree::groupOTU(tree, list_samples)
    p <- ggtree::ggtree(tree, ladderize=F) +
    	geom_tippoint(aes(color=group), size=size) +
    	scale_color_manual(values=label_colors, name=label)

    if(consensus){
    	p <- p + geom_tiplab(aes(color=group), size=5, hjust=-0.5)
    }

    p
  }

}

