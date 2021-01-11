#' identifies possible normal cells in the dataset based on coefficient of variation
#'
#' @param scCNA scCNA object
#' @param resolution Numeric. Threshold which will be used to detect normal cells. Lower values are more likely to be normal cells. Defaults to 0.05.
#'
#' @return Adds is_normal column to the scCNA metadata. Can be accessed with colData(scCNA)
#' @export
#'
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter bind_rows
#'
#' @examples
findNormalCells <- function(scCNA,
                            resolution = 0.05) {

  if (!is.numeric(resolution)) {
    stop("Resolution must be of class numeric")
  }

  `%!in%` <- Negate(`%in%`)

  if (is.null(SummarizedExperiment::colData(scCNA)$subclones)) {
    stop("findNormalCells needs clustering information. Please run findClusters().")
  }

  seg <- segment_ratios(scCNA) %>% t() %>% as.data.frame()

  seg_cl <- split(seg, as.character(SummarizedExperiment::colData(scCNA)$subclones))

  seg_cl_complete <- dplyr::bind_rows(seg_cl, .id = "cluster")

  cv <- lapply(seg_cl, function(x) apply(x, 1, function(z)  sd(z)/mean(z)))

  cv_cl_mean <- lapply(cv, mean)

  cv_df <- tibble::enframe(unlist(cv_cl_mean),
                           name = "cluster",
                           value = "CV")

  cv_df_low_cv <- cv_df %>%
    dplyr::filter(CV < resolution)

  message(paste0("Copykit detected ", nrow(cv_df_low_cv), " that are possibly normal cells with a resolution of: ", resolution))
  print(cv_df_low_cv)

  info <- data.frame(sample = rownames(seg_cl_complete),
                     cluster = seg_cl_complete$cluster,
                     stringsAsFactors = FALSE)

  info <- info %>%
    dplyr::mutate(is_normal = dplyr::case_when(
      cluster %!in% cv_df_low_cv$cluster ~ FALSE,
      TRUE ~ TRUE
    ))

  # reordering info to add to metadata
  info <- info[match(SummarizedExperiment::colData(scCNA)$sample, info$sample),]

  SummarizedExperiment::colData(scCNA)$is_normal <- info$is_normal

  message("Done. Information was added to metadata column 'is_normal'.")

  return(scCNA)

}


