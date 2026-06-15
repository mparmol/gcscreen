#' Call hits by FDR threshold and dual-method convergence
#'
#' Identifies compounds that pass an FDR threshold in both ANOVA and LMM
#' (the dual-convergence criterion used in the primary pipeline), and
#' optionally filters by a maximum AUC_rel threshold (growth inhibitors only).
#'
#' @param anova_res Data frame returned by [run_anova_screen()], with at least
#'   columns `Chemical`, `FDR`, and `AUC_rel`.
#' @param lmm_res Data frame returned by [run_lmm_screen()], with at least
#'   columns `Chemical` and `FDR`.
#' @param fdr_threshold Scalar FDR cutoff applied to both methods. Default
#'   `0.05`.
#' @param auc_rel_max Upper bound for AUC_rel: only compounds with
#'   `AUC_rel <= auc_rel_max` are considered hits (growth inhibitors). Set to
#'   `Inf` to include compounds with increased growth. Default `0.99`.
#' @param exclude Character vector of compound names to exclude regardless of
#'   statistics (e.g., optical artefacts). Default `NULL`.
#'
#' @return A data frame of confirmed hits with columns from `anova_res` plus
#'   `FDR_lmm` (from `lmm_res`), sorted by ascending `AUC_rel`.
#' @export
#' @examples
#' # Small demo (not a real screen)
#' anova_res <- data.frame(
#'   Chemical = c("DrugA", "DrugB", "DrugC", "DrugD"),
#'   Plate    = 1L,
#'   n_reps   = 3L,
#'   AUC_rel  = c(0.1, 0.5, 0.95, 0.2),
#'   p_anova  = c(0.001, 0.003, 0.6, 0.002),
#'   FDR      = c(0.004, 0.006, 0.6, 0.005),
#'   stringsAsFactors = FALSE
#' )
#' lmm_res <- data.frame(
#'   Chemical = c("DrugA", "DrugB", "DrugC", "DrugD"),
#'   FDR      = c(0.003, 0.04, 0.7, 0.06),
#'   stringsAsFactors = FALSE
#' )
#' call_hits(anova_res, lmm_res)
call_hits <- function(anova_res, lmm_res, fdr_threshold = 0.05,
                      auc_rel_max = 0.99, exclude = NULL) {
  # Compounds passing ANOVA FDR
  pass_anova <- anova_res$Chemical[
    !is.na(anova_res$FDR) & anova_res$FDR < fdr_threshold &
      !is.na(anova_res$AUC_rel) & anova_res$AUC_rel <= auc_rel_max
  ]

  # Compounds passing LMM FDR
  pass_lmm <- lmm_res$Chemical[
    !is.na(lmm_res$FDR) & lmm_res$FDR < fdr_threshold
  ]

  # Intersection (dual convergence)
  hits <- intersect(pass_anova, pass_lmm)

  # Remove excluded compounds
  if (!is.null(exclude)) hits <- setdiff(hits, exclude)

  if (length(hits) == 0L) {
    message("No hits identified at FDR < ", fdr_threshold)
    return(data.frame())
  }

  out <- anova_res[anova_res$Chemical %in% hits, ]
  lmm_fdr <- lmm_res[lmm_res$Chemical %in% hits, c("Chemical", "FDR")]
  names(lmm_fdr)[2L] <- "FDR_lmm"
  out <- merge(out, lmm_fdr, by = "Chemical", all.x = TRUE)
  out[order(out$AUC_rel), ]
}
