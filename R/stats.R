#' Per-compound ANOVA against plate controls
#'
#' Runs a one-way ANOVA comparing the replicate AUC values of one compound
#' against the negative-control wells of the same plate(s). Returns the raw
#' p-value.
#'
#' @param auc_compound Numeric vector of compound AUC values (one per
#'   biological replicate, typically 2–4).
#' @param auc_control  Numeric vector of control AUC values from the same
#'   plate(s) (typically 16 per plate-run × n_reps).
#' @param compound_name Character string used in error messages. Default `""`.
#'
#' @return Scalar p-value from the F-test, or `NA_real_` if the model could
#'   not be fitted (fewer than 2 observations, zero variance, etc.).
#' @export
#' @examples
#' ctrl <- rnorm(48, mean = 42, sd = 3)
#' cmpd <- rnorm(3, mean = 20, sd = 2)
#' screen_anova(cmpd, ctrl)
screen_anova <- function(auc_compound, auc_control, compound_name = "") {
  auc_compound <- as.numeric(auc_compound)
  auc_control  <- as.numeric(auc_control)
  n_cmpd <- sum(!is.na(auc_compound))
  n_ctrl <- sum(!is.na(auc_control))
  if (n_cmpd < 2L || n_ctrl < 2L) return(NA_real_)

  dat <- data.frame(
    auc   = c(auc_compound, auc_control),
    group = factor(c(rep("compound", n_cmpd), rep("control", n_ctrl)))
  )
  tryCatch({
    mod <- stats::lm(auc ~ group, data = dat)
    an  <- stats::anova(mod)
    an[["Pr(>F)"]][1L]
  }, error = function(e) {
    warning("ANOVA failed for '", compound_name, "': ", conditionMessage(e))
    NA_real_
  })
}

#' Per-compound linear mixed model with plate as random effect
#'
#' Fits `auc ~ Chemical + (1 | Plate)` using [lme4::lmer()] for compounds that
#' appear on two or more plates. Falls back to a standard ANOVA (equivalent to
#' [screen_anova()]) when fewer than two plates are available.
#'
#' @param data A data frame containing at least the columns specified by
#'   `auc_col`, `group_col`, and `plate_col`. Must include both the compound
#'   and control rows (labeled in `group_col`).
#' @param auc_col   Name of the column containing AUC values. Default `"auc"`.
#' @param group_col Name of the column identifying compound vs. control.
#'   Default `"Chemical"`. The control level must be `"CONTROL"`.
#' @param plate_col Name of the column containing the plate identifier.
#'   Default `"plate"`.
#' @param compound_name Character string for error messages. Default `""`.
#'
#' @return A named list with elements:
#'   \item{p_value}{Raw p-value for the compound effect.}
#'   \item{estimate}{Estimated AUC difference (compound – control).}
#'   \item{se}{Standard error of the estimate.}
#'   \item{method}{`"LMM"` or `"ANOVA"` depending on which model was used.}
#' @importFrom lmerTest lmer
#' @export
#' @examples
#' # Build a minimal long-format data frame
#' set.seed(42)
#' df <- data.frame(
#'   auc      = c(rnorm(3, 20, 2), rnorm(48, 42, 3)),
#'   Chemical = c(rep("DrugA", 3), rep("CONTROL", 48)),
#'   plate    = c(c(2, 3, 5), sample(2:6, 48, replace = TRUE))
#' )
#' screen_lmm(df, compound_name = "DrugA")
screen_lmm <- function(data, auc_col = "auc", group_col = "Chemical",
                        plate_col = "plate", compound_name = "") {
  data[[plate_col]] <- factor(data[[plate_col]])
  data[[group_col]] <- factor(data[[group_col]],
                               levels = c("CONTROL",
                                          setdiff(unique(data[[group_col]]),
                                                  "CONTROL")))

  n_plates <- length(unique(data[[plate_col]][data[[group_col]] != "CONTROL"]))

  na_result <- list(p_value = NA_real_, estimate = NA_real_,
                    se = NA_real_, method = NA_character_)

  if (nrow(data) < 4L) return(na_result)

  formula_lmm <- stats::as.formula(
    paste(auc_col, "~ Chemical + (1 |", plate_col, ")")
  )
  formula_lm <- stats::as.formula(
    paste(auc_col, "~", group_col)
  )

  if (n_plates < 2L) {
    # Fallback to ANOVA
    tryCatch({
      mod <- stats::lm(formula_lm, data = data)
      an  <- stats::anova(mod)
      cf  <- stats::coef(summary(mod))
      list(p_value  = an[["Pr(>F)"]][1L],
           estimate = if (nrow(cf) >= 2L) cf[2L, "Estimate"] else NA_real_,
           se       = if (nrow(cf) >= 2L) cf[2L, "Std. Error"] else NA_real_,
           method   = "ANOVA")
    }, error = function(e) {
      warning("ANOVA fallback failed for '", compound_name, "': ",
              conditionMessage(e))
      na_result
    })
  } else {
    tryCatch({
      fit <- lmerTest::lmer(formula_lmm, data = data, REML = FALSE)
      cf  <- stats::coef(summary(fit))
      if (nrow(cf) < 2L) return(na_result)
      list(p_value  = cf[2L, "Pr(>|t|)"],
           estimate = cf[2L, "Estimate"],
           se       = cf[2L, "Std. Error"],
           method   = "LMM")
    }, error = function(e) {
      # LMM singular or failed — fall back to ANOVA
      tryCatch({
        mod <- stats::lm(formula_lm, data = data)
        an  <- stats::anova(mod)
        cf  <- stats::coef(summary(mod))
        list(p_value  = an[["Pr(>F)"]][1L],
             estimate = if (nrow(cf) >= 2L) cf[2L, "Estimate"] else NA_real_,
             se       = if (nrow(cf) >= 2L) cf[2L, "Std. Error"] else NA_real_,
             method   = "ANOVA_fallback")
      }, error = function(e2) {
        warning("LMM and ANOVA both failed for '", compound_name, "'")
        na_result
      })
    })
  }
}

#' Run per-compound ANOVA for all compounds in a screen
#'
#' Applies [screen_anova()] to every compound in a tidy long-format data frame
#' and returns a results table with raw p-values and Benjamini-Hochberg FDR.
#'
#' @param screen_data Data frame with at least columns `Chemical`, `auc`,
#'   `plate`, and `run`. Control wells must be labeled `"CONTROL"` in the
#'   `Chemical` column.
#' @param fdr_method Correction method passed to [stats::p.adjust()].
#'   Default `"BH"`.
#'
#' @return A data frame with one row per compound and columns `Chemical`,
#'   `Plate`, `n_reps`, `AUC_rel` (mean relative AUC), `p_anova`, `FDR`.
#' @export
run_anova_screen <- function(screen_data, fdr_method = "BH") {
  compounds <- unique(screen_data$Chemical[screen_data$Chemical != "CONTROL"])
  ctrl_all  <- screen_data[screen_data$Chemical == "CONTROL", ]

  results <- lapply(compounds, function(cmpd) {
    sub_cmpd <- screen_data[screen_data$Chemical == cmpd, ]
    pl_id    <- sub_cmpd$plate[1L]
    sub_ctrl <- ctrl_all[ctrl_all$plate == pl_id, ]

    pval <- screen_anova(sub_cmpd$auc, sub_ctrl$auc, compound_name = cmpd)

    data.frame(
      Chemical  = cmpd,
      Plate     = pl_id,
      n_reps    = nrow(sub_cmpd),
      AUC_rel   = mean(sub_cmpd$AUC_relative, na.rm = TRUE),
      p_anova   = pval,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  out$FDR <- stats::p.adjust(out$p_anova, method = fdr_method)
  out[order(out$FDR), ]
}

#' Run per-compound LMM for all compounds in a screen
#'
#' Applies [screen_lmm()] to every compound and returns a results table with
#' raw p-values, estimates, and Benjamini-Hochberg FDR.
#'
#' @inheritParams run_anova_screen
#' @return A data frame with one row per compound and columns `Chemical`,
#'   `Plate`, `AUC_rel`, `Estimate`, `SE`, `p_lmm`, `FDR`, `method`.
#' @export
run_lmm_screen <- function(screen_data, fdr_method = "BH") {
  compounds <- unique(screen_data$Chemical[screen_data$Chemical != "CONTROL"])
  ctrl_all  <- screen_data[screen_data$Chemical == "CONTROL", ]

  results <- lapply(compounds, function(cmpd) {
    sub_cmpd <- screen_data[screen_data$Chemical == cmpd, ]
    pl_ids   <- unique(sub_cmpd$plate)
    sub_ctrl <- ctrl_all[ctrl_all$plate %in% pl_ids, ]

    long <- rbind(
      data.frame(Chemical = cmpd, auc = sub_cmpd$auc,
                 plate = sub_cmpd$plate, stringsAsFactors = FALSE),
      data.frame(Chemical = "CONTROL", auc = sub_ctrl$auc,
                 plate = sub_ctrl$plate, stringsAsFactors = FALSE)
    )

    res <- screen_lmm(long, compound_name = cmpd)

    data.frame(
      Chemical = cmpd,
      Plate    = pl_ids[1L],
      AUC_rel  = mean(sub_cmpd$AUC_relative, na.rm = TRUE),
      Estimate = res$estimate,
      SE       = res$se,
      p_lmm    = res$p_value,
      method   = res$method,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, results)
  out$FDR <- stats::p.adjust(out$p_lmm, method = fdr_method)
  out[order(out$FDR), ]
}
