#' Example 96-well kinetic screen dataset
#'
#' A simulated phenotypic screen with 3 plates (plate IDs 2, 3, 5) × 3
#' biological replicates, totalling 9 plate-runs. Each plate-run contains
#' 80 compound wells (columns 2–11) and 16 control wells (columns 1 and 12),
#' measured at 16 time points (0–60 h, 4 h resolution).
#'
#' Twelve "true hits" (4 per plate) were planted with carrying capacity
#' k ∈ (0.05, 0.35) to simulate growth inhibition. The remaining 228 compounds
#' have k ∈ (0.65, 0.95) (normal growth). A mild row-gradient positional
#' effect is included to exercise the B-score correction.
#'
#' @format A named list of 9 data frames (one per plate-run), each with columns:
#' \describe{
#'   \item{`run`}{Character. Plate-run identifier, e.g. `"placa_2_rep_A"`.}
#'   \item{`plate`}{Integer. Plate number (2, 3, or 5).}
#'   \item{`well_pos`}{Integer. Well position 1–96 (row-major).}
#'   \item{`time`}{Numeric. Time in hours (0, 4, 8, …, 60).}
#'   \item{`od`}{Numeric. Simulated OD600 measurement.}
#'   \item{`Chemical`}{Character. Compound name or `"CONTROL"`.}
#' }
#' @examples
#' data(example_screen)
#' names(example_screen)
#' head(example_screen[["placa_2_rep_A"]])
"example_screen"
