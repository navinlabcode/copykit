#' findAneuploidCells
#'
#' Find cells that are not aneuploid in the dataset.
#'
#' @param scCNA The CopyKit object
#' @param assay String with the name of the assay to pull data from to find
#' normal cells.
#' @param resolution A numeric scalar used as threshold to detect normal cells.
#' @param remove_XY A boolean that removes chrX and chrY from the analysis.
#' @param simul A boolean that if TRUE adds a simulated normal dataset to boost
#' identifying normal cells in datasets with small proportions of normal cells.
#' @param seed Seed passed on to reproduce simulated CV of normal cells.
#'
#' @details performs a sample-wise calculation of the segment means coefficient
#'  of variation and fits a Gaussian mixture model to the observed distribution f
#'  rom all cells. To increase the sensitivity of the model, the expected
#'  distribution of the coefficient of variation for euploid cells is simulated
#'  for a thousand cells (mean = 0, sd = 0.01). This way, CopyKit can adequately
#'  detect euploid cells even in datasets with limited amounts of euploid cells
#'  and guarantees that no aneuploid cell will be removed from datasets without
#'  any euploid cells. The distribution with the smallest coefficient of variance
#'  is assumed originate from normal cells. Cells are classified as euploid
#'  if they have a coefficient of variance smaller than the mean plus five times
#'  the standard deviation of the normal cell distribution.
#'
#' @return information is added to \code{\link[SummarizedExperiment]{colData}}
#' in a columns named 'is_aneuploid' being TRUE if a cell is detected as
#' aneuploid and FALSE if the cell is detected as euploid.
#'
#' @export
#'
#' @importFrom tibble enframe
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr filter bind_rows case_when
#' @importFrom mixtools normalmixEM
#' @importFrom stats rnorm sd
#'
#' @examples
#' copykit_obj <- copykit_example()
#' copykit_obj <- findAneuploidCells(copykit_obj)
findAneuploidCells <- function(scCNA,
    assay = "segment_ratios",
    resolution = "auto",
    remove_XY = TRUE,
    simul = TRUE,
    seed = 17) {

    # bindings for NSE (non-standard evaluation)
    is_aneuploid <- NULL

    if (remove_XY == FALSE & simul == TRUE) {
        stop("Argument simul can't be used if remove_XY == FALSE.")
    }

    if (resolution != "auto" & !is.numeric(resolution)) {
        stop("Resolution must be of class numeric")
    }

    # retrieving data
    rg <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
    seg <- SummarizedExperiment::assay(scCNA, assay)
    ncells <- ncol(scCNA)

    if (remove_XY == TRUE) {
        rg <- rg %>%
            dplyr::filter(
                !stringr::str_detect(seqnames, "X"),
                !stringr::str_detect(seqnames, "Y")
            )

        seg <- seg[1:nrow(rg), ]
    }

    # calculating the coefficient of variation
    cv <- vapply(
        seg, function(z) {
            sd(z) / mean(z)
        },
        numeric(1)
    )

    if (simul == TRUE) {
        withr::with_seed(seed,
            cv_simul <- rnorm(ncells,
                              mean = 0,
                              sd = 0.01
            )
        )

        names(cv_simul) <- paste0("simul", 1:length(cv_simul))

        cv <- c(cv_simul, cv)
    }

    if (resolution == "auto") {
        fit <- tryCatch(
            mixtools::normalmixEM(cv),
            error = function(e) {
                message("Could not identify aneuploid cells in the dataset.")
                message("Marking all cells as diploid.")
                message("Check colData(scCNA)$find_normal_cv.")
                return("error")
            }
        )

        # determining resolution
        if (length(fit) > 1) {
            resolution <- fit$mu[1] + 5 * fit$sigma[1]
        } else {
            resolution <- 1
        }
    }

    if (simul == TRUE) {
        cv <- cv[!grepl("simul", names(cv))]
    }

    cv_df <- tibble::enframe(cv,
        name = "sample",
        value = "CV"
    )

    cv_df_low_cv <- cv_df %>%
        dplyr::mutate(is_aneuploid = case_when(
            CV > resolution ~ TRUE,
            TRUE ~ FALSE
        ))

    message(
        "Copykit detected ",
        nrow(cv_df_low_cv %>%
            dplyr::filter(is_aneuploid == FALSE)),
        " that are possibly euploid cells using a resolution of: ",
        round(resolution, 3)
    )

    # reordering info to add to metadata
    info <-
        cv_df_low_cv[match(
            SummarizedExperiment::colData(scCNA)$sample,
            cv_df_low_cv$sample
        ), ]

    SummarizedExperiment::colData(scCNA)$is_aneuploid <- info$is_aneuploid
    SummarizedExperiment::colData(scCNA)$find_normal_cv <-
        round(info$CV, 2)

    message("Added information to colData(CopyKit).")

    return(scCNA)
}
