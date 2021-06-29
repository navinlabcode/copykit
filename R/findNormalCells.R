#' findNormalCells
#'
#' Find cells that are not aneuploid in the dataset.
#'
#' @param scCNA scCNA object
#' @param assay String with the name of the assay to pull data from to find normal cells.
#' @param resolution A numeric scalar used as threshold to detect normal cells. See details.
#' @param remove_XY A boolean that removes chrX and chrY from the analysis. Recommended.
#' @param simul A boolean that if TRUE adds a simulated normal dataset to boost
#' identifying normal cells in datasets with small proportions of normal cells.
#'
#' @details performs a sample-wise calculation of the segment means coefficient
#'  of variation and fits a normal mixture model to the observed distribution f
#'  rom all cells. To increase the sensitivity of the model, the expected
#'  distribution of the coefficient of variation for diploid cells is simulated
#'  for a thousand cells (mean = 0, sd = 0.01). This way, CopyKit can adequately
#'  detect normal cells even in datasets with limited amounts of diploid cells
#'  and guarantees that no aneuploid cell will be removed from datasets without
#'  any normal cells. The distribution with the smallest coefficient of variance
#'  is assumed to be originating from normal cells. Cells are classified as normal
#'  if they have a coefficient of variance smaller than the mean plus five times
#'  the standard deviation of the normal cell distribution.
#'
#' @return information is added to \code{\link[SummarizedExperiment]{colData}}
#' in a columns named 'is_normal' being TRUE if a cell is detected as normal and
#' FALSE if the cell is detected as aneuploid.
#'
#' @export
#'
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter bind_rows case_when
#' @importFrom mixtools normalmixEM
#'
#' @examples
findNormalCells <- function(scCNA,
                            assay = 'segment_ratios',
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
  seg <- SummarizedExperiment::assay(scCNA, assay)

  if (remove_XY == TRUE) {
    rg <- rg %>%
      dplyr::filter(!stringr::str_detect(seqnames, "X"),
                    !stringr::str_detect(seqnames, "Y"))

    seg <- seg[1:nrow(rg), ]
  }

  # calculating the coeficient of variation
  cv <-  sapply(seg, function(z)
    sd(z) / mean(z))

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
    resolution <- fit$mu[1] + 5 * fit$sigma[1]
  }

  if (simul == TRUE) {
    cv <- cv[!grepl("simul", names(cv))]
  }

  cv_df <- tibble::enframe(cv,
                           name = "sample",
                           value = "CV")

  cv_df_low_cv <- cv_df %>%
    dplyr::mutate(is_normal = case_when(CV > resolution ~ FALSE,
                                        TRUE ~ TRUE))

  message(
    paste0(
      "Copykit detected ",
      nrow(cv_df_low_cv %>%
             dplyr::filter(is_normal == TRUE)),
      " that are possibly normal cells using a resolution of: ",
      round(resolution, 3)
    )
  )

  # reordering info to add to metadata
  info <-
    cv_df_low_cv[match(SummarizedExperiment::colData(scCNA)$sample,
                       cv_df_low_cv$sample), ]

  SummarizedExperiment::colData(scCNA)$is_normal <- info$is_normal
  SummarizedExperiment::colData(scCNA)$find_normal_cv <-
    round(info$CV, 2)

  message("Done. Information was added to metadata column 'is_normal'.")

  return(scCNA)

}
