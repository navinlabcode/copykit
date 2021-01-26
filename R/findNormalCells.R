#' identifies possible normal cells in the dataset based on coefficient of variation
#'
#' @param scCNA scCNA object
#' @param resolution Numeric. Threshold which will be used to detect normal cells.
#'
#' @return Adds is_normal column to the scCNA metadata. Can be accessed with colData(scCNA)
#' @export
#'
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter bind_rows
#' @importFrom mixtools normalmixEM
#'
#' @examples
findNormalCells <- function(scCNA,
                            resolution = "fit") {

  if (resolution != "fit" & !is.numeric(resolution)) {
    stop("Resolution must be of class numeric")
  }

  seg <- segment_ratios(scCNA)

  cv <-  sapply(seg, function(z) sd(z)/mean(z))

  if (resolution == "fit") {
    fit <- mixtools::normalmixEM(cv)
    resolution <- fit$mu[1] + 5*fit$sigma[1]
  }

  cv_df <- tibble::enframe(cv,
                           name = "sample",
                           value = "CV")

  cv_df_low_cv <- cv_df %>%
    dplyr::mutate(is_normal = case_when(
      CV > resolution ~ FALSE,
      TRUE ~ TRUE)
    )

  message(paste0("Copykit detected ",
                 nrow(cv_df_low_cv %>%
                        dplyr::filter(is_normal == TRUE)),
                 " that are possibly normal cells using a resolution of: ",
                 resolution))

  # reordering info to add to metadata
  info <- cv_df_low_cv[match(SummarizedExperiment::colData(scCNA)$sample, cv_df_low_cv$sample),]

  SummarizedExperiment::colData(scCNA)$is_normal <- info$is_normal
  SummarizedExperiment::colData(scCNA)$find_normal_cv <- round(info$CV,2)

  message("Done. Information was added to metadata column 'is_normal'.")

  return(scCNA)

}
