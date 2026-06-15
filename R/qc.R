#' Coefficient of variation of control wells
#'
#' @param ctrl_auc Numeric vector of control well AUC values.
#' @param na.rm Logical; remove `NA` values. Default `TRUE`.
#' @return CV as a percentage (SD / mean × 100).
#' @export
#' @examples
#' ctrl_cv(c(40, 42, 41, 43, 44, 40, 39, 42, 41, 43, 44, 42, 40, 41, 43, 42))
ctrl_cv <- function(ctrl_auc, na.rm = TRUE) {
  m <- mean(ctrl_auc, na.rm = na.rm)
  s <- stats::sd(ctrl_auc, na.rm = na.rm)
  if (is.na(m) || m == 0) return(NA_real_)
  100 * s / m
}

#' Z'-factor
#'
#' Computes the Z'-factor (Zhang et al., 1999) to assess assay quality.
#' Z' = 1 - 3*(sd_pos + sd_neg) / |mean_pos - mean_neg|.
#'
#' @param pos_auc Numeric vector of positive-control (or hit) AUC values.
#' @param neg_auc Numeric vector of negative-control AUC values.
#' @param na.rm Logical. Default `TRUE`.
#' @return Scalar Z'-factor. Returns `NA` if fewer than 2 observations in
#'   either group or if the mean difference is zero.
#' @references Zhang J-H, Chung TDY, Oldenburg KR (1999). A Simple Statistical
#'   Parameter for Use in Evaluation and Validation of High Throughput Screening
#'   Assays. *J Biomol Screen* 4(2):67-73.
#' @export
zprime <- function(pos_auc, neg_auc, na.rm = TRUE) {
  pos_auc <- stats::na.omit(pos_auc)
  neg_auc <- stats::na.omit(neg_auc)
  if (length(pos_auc) < 2L || length(neg_auc) < 2L) return(NA_real_)
  denom <- abs(mean(pos_auc) - mean(neg_auc))
  if (denom == 0) return(NA_real_)
  1 - 3 * (stats::sd(pos_auc) + stats::sd(neg_auc)) / denom
}

#' Robust Z'-factor
#'
#' A robust version of the Z'-factor using medians and MADs instead of means
#' and SDs. RZ' = 1 - 3*(MAD_pos + MAD_neg) / |median_pos - median_neg|.
#'
#' @inheritParams zprime
#' @return Scalar RZ'-factor.
#' @export
rz_prime <- function(pos_auc, neg_auc, na.rm = TRUE) {
  pos_auc <- stats::na.omit(pos_auc)
  neg_auc <- stats::na.omit(neg_auc)
  if (length(pos_auc) < 2L || length(neg_auc) < 2L) return(NA_real_)
  denom <- abs(stats::median(pos_auc) - stats::median(neg_auc))
  if (denom == 0) return(NA_real_)
  1 - 3 * (stats::mad(pos_auc) + stats::mad(neg_auc)) / denom
}

#' Strictly Standardized Mean Difference (SSMD)
#'
#' SSMD = (mean_pos - mean_neg) / sqrt(sd_pos^2 + sd_neg^2).
#' A value ≥ 3 indicates excellent assay quality (Zhang, 2007).
#'
#' @inheritParams zprime
#' @return Scalar SSMD.
#' @references Zhang XD (2007). A pair of new statistical parameters for
#'   quality control in RNA interference high-throughput screening assays.
#'   *Genomics* 89(4):552-61.
#' @export
ssmd <- function(pos_auc, neg_auc, na.rm = TRUE) {
  pos_auc <- stats::na.omit(pos_auc)
  neg_auc <- stats::na.omit(neg_auc)
  if (length(pos_auc) < 2L || length(neg_auc) < 2L) return(NA_real_)
  denom <- sqrt(stats::var(pos_auc) + stats::var(neg_auc))
  if (is.na(denom) || denom == 0) return(NA_real_)
  (mean(pos_auc) - mean(neg_auc)) / denom
}

#' Plate-run quality metrics
#'
#' Computes CV, Z'-factor, robust Z'-factor (RZ'), and SSMD for a single
#' plate-run, given control and (optionally) hit AUC values.
#'
#' @param ctrl_auc   Numeric vector of negative-control AUC values.
#' @param hit_auc    Numeric vector of hit (positive-control surrogate) AUC
#'   values. If `NULL` or fewer than 2 observations, Z', RZ', and SSMD
#'   are returned as `NA`.
#' @param cv_warn    CV threshold (%) above which a warning flag is set.
#'   Default `20`.
#'
#' @return A named list with elements `cv`, `zprime`, `rz_prime`, `ssmd`,
#'   and `cv_flag` (`"OK"` or `"HIGH_CV"`).
#' @export
#' @examples
#' ctrl <- rnorm(16, mean = 42, sd = 3)
#' hits <- rnorm(5, mean = 10, sd = 2)
#' plate_qc(ctrl, hits)
plate_qc <- function(ctrl_auc, hit_auc = NULL, cv_warn = 20) {
  cv_val <- ctrl_cv(ctrl_auc)
  if (!is.null(hit_auc) && sum(!is.na(hit_auc)) >= 2L) {
    zp  <- zprime(hit_auc, ctrl_auc)
    rzp <- rz_prime(hit_auc, ctrl_auc)
    sm  <- ssmd(hit_auc, ctrl_auc)
  } else {
    zp <- rzp <- sm <- NA_real_
  }
  list(
    cv       = cv_val,
    zprime   = zp,
    rz_prime = rzp,
    ssmd     = sm,
    cv_flag  = ifelse(!is.na(cv_val) && cv_val > cv_warn, "HIGH_CV", "OK")
  )
}
