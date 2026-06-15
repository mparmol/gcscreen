#' Convert well position to row and column indices
#'
#' Maps an integer well position (1–96) in a standard 8×12 microplate layout
#' (row-major order: well 1 = A1, well 12 = A12, well 13 = B1, ..., 96 = H12)
#' to numeric row (1–8) and column (1–12) indices.
#'
#' @param well_pos Integer vector of well positions (1–96).
#' @return A data frame with columns `row_idx` (1–8) and `col_idx` (1–12).
#' @export
#' @examples
#' well_to_rowcol(c(1, 12, 13, 96))
well_to_rowcol <- function(well_pos) {
  well_pos <- as.integer(well_pos)
  if (any(well_pos < 1 | well_pos > 96, na.rm = TRUE))
    stop("well_pos must be between 1 and 96")
  data.frame(
    row_idx = ceiling(well_pos / 12L),
    col_idx = ((well_pos - 1L) %% 12L) + 1L
  )
}

#' Trapezoidal area under the curve
#'
#' Computes the area under a curve using the trapezoidal rule.
#'
#' @param time Numeric vector of time points (must be strictly increasing).
#' @param od   Numeric vector of OD600 values, same length as `time`.
#' @return Scalar: trapezoidal AUC.
#' @export
#' @examples
#' t <- seq(0, 60, by = 4)
#' od <- plogis((t - 20) / 5) * 0.8
#' auc_trapz(t, od)
auc_trapz <- function(time, od) {
  n <- length(time)
  if (n != length(od)) stop("'time' and 'od' must have the same length")
  if (n < 2L) return(NA_real_)
  sum(diff(time) * (od[-n] + od[-1L])) / 2
}

#' Baseline-correct OD600 time series
#'
#' Subtracts the value at t = 0 (first time point) from every subsequent
#' measurement in each well, then floors negative values at zero.
#'
#' @param od_matrix Numeric matrix where rows are time points and columns are
#'   wells. The first row must correspond to t = 0.
#' @return Matrix of the same dimensions with baseline-corrected, non-negative
#'   values.
#' @export
baseline_correct <- function(od_matrix) {
  if (!is.matrix(od_matrix) && !is.data.frame(od_matrix))
    stop("'od_matrix' must be a matrix or data frame")
  od_matrix <- as.matrix(od_matrix)
  baseline  <- od_matrix[1L, , drop = FALSE]
  corrected <- sweep(od_matrix, 2L, baseline, "-")
  corrected[corrected < 0] <- 0
  corrected
}
