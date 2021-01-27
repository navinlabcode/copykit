#' identifies possible normal cells in the dataset based on coefficient of variation
#'
#' @param scCNA scCNA object
#' @param resolution Numeric. Threshold which will be used to detect normal cells.
#' @param remove_XY Boolean. Removes chrX and chrY from the analysis. Recommended.
#' @param simul Add a simulated normal dataset to boost identifying normal cells when a dataset has a small proportion of those.
#'
#' @return Adds is_normal column to the scCNA metadata. Can be accessed with colData(scCNA)
#' @export
#'
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter bind_rows case_when
#' @importFrom mixtools normalmixEM
#'
#' @examples
findNormalCells <- function(scCNA,
                            resolution = "auto",
                            remove_XY = TRUE,
                            simul = TRUE) {

  if (remove_XY == FALSE & simul == TRUE) {
    stop("Argument simul can't be used if remove_XY == FALSE.")
  }

  if (resolution != "auto" & !is.numeric(resolution)) {
    stop("Resolution must be of class numeric")
  }

  # retrieving data
  rg <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
  seg <- segment_ratios(scCNA)

  if (remove_XY == TRUE) {
    rg <- rg %>%
      dplyr::filter(!stringr::str_detect(seqnames, "X"),
                    !stringr::str_detect(seqnames, "Y"))

    seg <- seg[1:nrow(rg),]
  }

  # calculating the coeficient of variation
  cv <-  sapply(seg, function(z) sd(z)/mean(z))

  if (simul == TRUE) {
    set.seed(17)
    cv_simul <- rnorm(1000,
          mean = 0,
          sd = 0.01)
    names(cv_simul) <- paste0("simul", 1:length(cv_simul))

    cv <- c(cv_simul, cv)

  }

  if (resolution == "auto") {
    fit <- mixtools::normalmixEM(cv)
    resolution <- fit$mu[1] + 5*fit$sigma[1]
  }

  if (simul == TRUE) {
    cv <- cv[!grepl("simul", names(cv))]
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
                 round(resolution,3)))

  # reordering info to add to metadata
  info <- cv_df_low_cv[match(SummarizedExperiment::colData(scCNA)$sample, cv_df_low_cv$sample),]

  SummarizedExperiment::colData(scCNA)$is_normal <- info$is_normal
  SummarizedExperiment::colData(scCNA)$find_normal_cv <- round(info$CV,2)

  message("Done. Information was added to metadata column 'is_normal'.")

  return(scCNA)

}
