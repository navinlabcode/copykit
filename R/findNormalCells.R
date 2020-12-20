#' identifies possible normal cells in the dataset based on coefficient of variation
#'
#' @param scCNA scCNA object
#'
#' @return Adds is_normal column to the scCNA metadata. Can be accessed with colData(scCNA)
#' @export
#'
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter
#'
#' @examples
findNormalCells <- function(scCNA) {

  if (is.null(SummarizedExperiment::colData(scCNA)$subclones)) {
    stop("findNormalCells needs clustering information. Please run findClusters().")
  }

  seg <- segment_ratios(scCNA) %>% t() %>% as.data.frame()

  seg_cl <- split(seg, as.character(colData(scCNA)$subclones))

  seg_cl_complete <- bind_rows(seg_cl, .id = "cluster")

  cv <- lapply(seg_cl, function(x) apply(x, 1, function(z)  sd(z)/mean(z)))

  cv_cl_mean <- lapply(cv, mean)

  cv_df <- tibble::enframe(unlist(cv_cl_mean),
                           name = "cluster",
                           value = "CV")

  cv_df <- cv_df %>%
    dplyr::filter(CV == min(CV))

  info <- data.frame(sample = rownames(seg_cl_complete),
                     cluster = seg_cl_complete$cluster,
                     stringsAsFactors = FALSE)

  info <- info %>%
    dplyr::mutate(is_normal = dplyr::case_when(
      cluster != cv_df$cluster ~ FALSE,
      TRUE ~ TRUE
    ))

  # reordering info to add to metadata
  info <- info[match(SummarizedExperiment::colData(scCNA)$sample, info$sample),]

  SummarizedExperiment::colData(scCNA)$is_normal <- info$is_normal

  return(scCNA)

}


