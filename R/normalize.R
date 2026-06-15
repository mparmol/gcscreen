#' Normalize compound AUC values to the control median
#'
#' Divides each compound's AUC by the median AUC of the negative-control wells
#' in the same plate-run, producing a dimensionless relative AUC (AUC_rel) where
#' 1.0 represents unimpaired growth.
#'
#' @param auc_compound Numeric vector of AUC values for compound wells in one
#'   plate-run.
#' @param auc_control  Numeric vector of AUC values for control wells in the
#'   same plate-run.
#' @param na.rm Logical; should `NA` values be removed when computing the
#'   control median? Default `TRUE`.
#'
#' @return Numeric vector of relative AUC values, same length as
#'   `auc_compound`. Returns `NA` if the control median is zero or `NA`.
#' @export
#' @examples
#' ctrl <- c(42, 44, 41, 43, 40, 45, 42, 43, 41, 44, 42, 43,
#'           41, 44, 43, 42)
#' cmpd <- c(38, 20, 5, 42, 41, 39)
#' normalize_controls(cmpd, ctrl)
normalize_controls <- function(auc_compound, auc_control, na.rm = TRUE) {
  med <- stats::median(auc_control, na.rm = na.rm)
  if (is.na(med) || med == 0) {
    warning("Control median is zero or NA; returning NA for all compound values")
    return(rep(NA_real_, length(auc_compound)))
  }
  as.numeric(auc_compound) / med
}
