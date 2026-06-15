#' B-score positional-bias correction for one plate-run
#'
#' Applies median polish (B-score) to the 8×10 sub-matrix of compound wells
#' (columns 2–11 of a standard 96-well plate) to remove row and column
#' positional effects caused by evaporation gradients, optical artefacts, or
#' thermal heterogeneity. Control wells (columns 1 and 12) are **not** modified.
#'
#' The correction subtracts the estimated row effect and column effect from each
#' compound AUC value while preserving the plate-wide overall level:
#' `AUC_corrected = AUC_raw - row_effect - col_effect`.
#'
#' @param well_pos Integer vector of well positions (1–96) for every well in
#'   the plate-run.
#' @param auc Numeric vector of AUC values, same order as `well_pos`.
#' @param ctrl_cols Integer vector of column indices that contain control wells
#'   and should **not** be corrected. Default `c(1L, 12L)`.
#' @param max_iter Maximum iterations for `stats::medpolish`. Default 10.
#'
#' @return Numeric vector of corrected AUC values, same length and order as
#'   the input `auc`. Control-well values are returned unchanged.
#' @export
#' @examples
#' # Simulate one plate-run: 96 wells, compound cols 2-11, ctrl cols 1 & 12
#' set.seed(1)
#' n    <- 96
#' wpos <- 1:n
#' auc_raw <- rnorm(n, mean = 40, sd = 4)
#' # Add a synthetic row gradient
#' rc <- well_to_rowcol(wpos)
#' auc_raw <- auc_raw + rc$row_idx * 0.5
#' auc_corr <- bscore_plate(wpos, auc_raw)
#' cor(well_to_rowcol(wpos)$row_idx[rc$col_idx %in% 2:11],
#'     auc_corr[rc$col_idx %in% 2:11])  # should be near 0
bscore_plate <- function(well_pos, auc, ctrl_cols = c(1L, 12L),
                         max_iter = 10L) {
  well_pos <- as.integer(well_pos)
  auc      <- as.numeric(auc)
  n        <- length(well_pos)
  if (length(auc) != n)
    stop("'well_pos' and 'auc' must have the same length")

  rc       <- well_to_rowcol(well_pos)
  is_cmpd  <- !(rc$col_idx %in% ctrl_cols)
  cmpd_idx <- which(is_cmpd)

  if (length(cmpd_idx) < 4L) {
    warning("Too few compound wells for B-score; returning uncorrected values")
    return(auc)
  }

  n_rows <- 8L
  n_cols <- 10L
  mat <- matrix(NA_real_, nrow = n_rows, ncol = n_cols,
                dimnames = list(seq_len(n_rows),
                                setdiff(seq_len(12L), ctrl_cols)))

  for (k in cmpd_idx) {
    r <- rc$row_idx[k]
    c <- rc$col_idx[k] - min(setdiff(seq_len(12L), ctrl_cols)) + 1L
    if (r >= 1L && r <= n_rows && c >= 1L && c <= n_cols)
      mat[r, c] <- auc[k]
  }

  mp <- tryCatch(
    stats::medpolish(mat, na.rm = TRUE, trace.iter = FALSE,
                     maxiter = max_iter),
    error = function(e) {
      warning("medpolish failed: ", conditionMessage(e),
              "; returning uncorrected values")
      NULL
    }
  )

  if (is.null(mp)) return(auc)

  auc_out <- auc
  for (k in cmpd_idx) {
    r <- rc$row_idx[k]
    c <- rc$col_idx[k] - min(setdiff(seq_len(12L), ctrl_cols)) + 1L
    if (r >= 1L && r <= n_rows && c >= 1L && c <= n_cols) {
      row_eff <- if (!is.na(mp$row[r])) mp$row[r] else 0
      col_eff <- if (!is.na(mp$col[c])) mp$col[c] else 0
      auc_out[k] <- auc[k] - row_eff - col_eff
    }
  }
  auc_out
}
