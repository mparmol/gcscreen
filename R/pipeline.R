#' Run the complete gcscreen pipeline on a list of plate-runs
#'
#' High-level wrapper that applies the full pipeline to a named list of
#' plate-run data frames: (1) baseline correction, (2) AUC calculation,
#' (3) B-score positional correction, (4) control normalization, (5) QC
#' metrics, (6) ANOVA per compound, (7) LMM per compound, (8) FDR-based
#' dual-convergence hit calling.
#'
#' @param platerun_list A named list of data frames, one per plate-run. Each
#'   data frame must have:
#'   \describe{
#'     \item{`time`}{Numeric time column (hours).}
#'     \item{`well_pos`}{Integer well position (1–96).}
#'     \item{`od`}{Numeric OD600 measurement.}
#'     \item{`Chemical`}{Character compound name or `"CONTROL"`.}
#'     \item{`plate`}{Integer plate identifier.}
#'   }
#'   Names of the list are used as `run` identifiers.
#' @param time_max Maximum time (hours) to include. Measurements beyond this
#'   are discarded. Default `60`.
#' @param ctrl_cols Integer vector of column positions (1–12) containing
#'   control wells. Default `c(1L, 12L)`.
#' @param ctrl_label Character label for control wells in the `Chemical`
#'   column. Default `"CONTROL"`.
#' @param plates_exclude Integer vector of plate IDs to exclude (e.g., those
#'   with CV > 45%). Default `NULL`.
#' @param cv_warn CV threshold (%) for flagging plate-runs. Default `20`.
#' @param fdr_threshold FDR cutoff for both ANOVA and LMM. Default `0.05`.
#' @param auc_rel_max AUC_rel upper bound for calling a hit (growth
#'   inhibitors). Default `0.99`.
#' @param exclude Character vector of compound names to exclude from hits
#'   regardless of statistics. Default `NULL`.
#'
#' @return An invisible list with elements:
#'   \describe{
#'     \item{`screen_data`}{Tidy long-format data frame with one row per
#'       well × plate-run, including AUC, AUC_relative, and B-score
#'       corrected values.}
#'     \item{`qc`}{Data frame of QC metrics per plate-run (CV, Z', SSMD,
#'       RZ').}
#'     \item{`anova_results`}{Per-compound ANOVA results with FDR.}
#'     \item{`lmm_results`}{Per-compound LMM results with FDR.}
#'     \item{`hits`}{Data frame of confirmed hits (ANOVA ∩ LMM, FDR <
#'       `fdr_threshold`).}
#'   }
#' @export
#' @examples
#' # Use the built-in example dataset
#' data(example_screen)
#' res <- run_screen(example_screen)
#' res$hits
run_screen <- function(platerun_list,
                       time_max       = 60,
                       ctrl_cols      = c(1L, 12L),
                       ctrl_label     = "CONTROL",
                       plates_exclude = NULL,
                       cv_warn        = 20,
                       fdr_threshold  = 0.05,
                       auc_rel_max    = 0.99,
                       exclude        = NULL) {

  if (!is.list(platerun_list) || is.null(names(platerun_list)))
    stop("'platerun_list' must be a named list of data frames")

  # Exclude requested plates
  if (!is.null(plates_exclude)) {
    keep <- sapply(platerun_list, function(df) {
      !(df$plate[1L] %in% plates_exclude)
    })
    platerun_list <- platerun_list[keep]
    message(sprintf("Excluded %d plate-run(s): plates %s",
                    sum(!keep),
                    paste(plates_exclude, collapse = ", ")))
  }

  if (length(platerun_list) == 0L) stop("No plate-runs remain after exclusion")

  all_data <- NULL
  qc_rows  <- NULL

  for (run_id in names(platerun_list)) {
    df <- platerun_list[[run_id]]

    # Truncate to time_max
    df <- df[df$time <= time_max, ]

    # AUC per well (trapezoidal, after baseline correction per well)
    well_aucs <- tapply(seq_len(nrow(df)), df$well_pos, function(idx) {
      sub <- df[idx, ]
      sub <- sub[order(sub$time), ]
      od_corr <- sub$od - sub$od[1L]
      od_corr[od_corr < 0] <- 0
      auc_trapz(sub$time, od_corr)
    })

    auc_df <- data.frame(
      run      = run_id,
      plate    = df$plate[1L],
      well_pos = as.integer(names(well_aucs)),
      auc      = as.numeric(well_aucs),
      Chemical = vapply(as.integer(names(well_aucs)), function(wp) {
        rows <- df[df$well_pos == wp, ]
        rows$Chemical[1L]
      }, character(1L)),
      stringsAsFactors = FALSE
    )

    # B-score on compound wells
    auc_df$auc_corrected <- bscore_plate(auc_df$well_pos, auc_df$auc,
                                          ctrl_cols = ctrl_cols)

    # Control normalization
    ctrl_auc  <- auc_df$auc_corrected[auc_df$Chemical == ctrl_label]
    med_ctrl  <- stats::median(ctrl_auc, na.rm = TRUE)
    auc_df$AUC_relative <- auc_df$auc_corrected / med_ctrl

    # QC
    hit_names <- if (!is.null(exclude)) {
      unique(auc_df$Chemical[auc_df$Chemical != ctrl_label &
                               !auc_df$Chemical %in% exclude])
    } else {
      unique(auc_df$Chemical[auc_df$Chemical != ctrl_label])
    }
    # Use pre-existing hits if available (first pass: all compounds as proxies)
    qc_res <- plate_qc(ctrl_auc, cv_warn = cv_warn)
    qc_rows <- rbind(qc_rows, data.frame(
      run      = run_id,
      plate    = df$plate[1L],
      n_ctrl   = length(ctrl_auc),
      cv       = qc_res$cv,
      zprime   = qc_res$zprime,
      rz_prime = qc_res$rz_prime,
      ssmd     = qc_res$ssmd,
      cv_flag  = qc_res$cv_flag,
      stringsAsFactors = FALSE
    ))

    all_data <- rbind(all_data, auc_df)
  }

  message(sprintf(
    "Processed %d plate-run(s) | %d unique compounds",
    length(platerun_list),
    length(unique(all_data$Chemical[all_data$Chemical != ctrl_label]))
  ))

  # Statistical testing
  message("Running ANOVA per compound...")
  anova_res <- run_anova_screen(all_data, fdr_method = "BH")

  message("Running LMM per compound...")
  lmm_res <- run_lmm_screen(all_data, fdr_method = "BH")

  # Hit calling
  hits <- call_hits(anova_res, lmm_res,
                    fdr_threshold = fdr_threshold,
                    auc_rel_max   = auc_rel_max,
                    exclude       = exclude)

  message(sprintf("Hits identified: %d", nrow(hits)))

  invisible(list(
    screen_data   = all_data,
    qc            = qc_rows,
    anova_results = anova_res,
    lmm_results   = lmm_res,
    hits          = hits
  ))
}
