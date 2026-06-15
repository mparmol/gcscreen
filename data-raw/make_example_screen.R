## Run from package root: source("data-raw/make_example_screen.R")
setwd("d:/Projects/Proyecto_curvas/gcscreen")

well_col <- function(wp) ((as.integer(wp) - 1L) %% 12L) + 1L

set.seed(2026)
time_pts <- seq(0, 60, by = 4)
n_time   <- length(time_pts)

true_hits <- list(
  "2" = c("Cycloheximide", "Lanoconazole", "Thimerosal", "Nitroxoline"),
  "3" = c("Amiodarone", "Hexachlorophene", "Clioquinol", "Clomiphene"),
  "5" = c("Doxorubicin", "Dequalinium", "Oxiconazole", "Tolcapone")
)

logistic_od <- function(t, k, t_mid = 18, sigma = 0.005) {
  od <- k / (1 + ((k - 0.02) / 0.02) * exp(-0.12 * (t - t_mid)))
  pmax(od + rnorm(length(t), 0, sigma), 0)
}

ctrl_pos <- (1:96)[well_col(1:96) %in% c(1L, 12L)]  # 16 wells
cmpd_pos <- (1:96)[well_col(1:96) %in% 2:11]         # 80 wells

example_screen <- list()

for (pl in c(2L, 3L, 5L)) {
  hits_pl    <- true_hits[[as.character(pl)]]
  cmpd_names <- c(hits_pl,
                  paste0("Compound_", sprintf("%03d", seq_len(80L - length(hits_pl)))))

  for (rep_letter in c("A", "B", "C")) {
    run_id <- sprintf("placa_%d_rep_%s", pl, rep_letter)
    all_rows <- vector("list", length(ctrl_pos) + length(cmpd_pos))
    idx <- 1L

    for (wp in ctrl_pos) {
      od <- logistic_od(time_pts, k = 0.8 + rnorm(1, 0, 0.03),
                        t_mid = 18 + rnorm(1, 0, 1))
      all_rows[[idx]] <- data.frame(
        run      = run_id,
        plate    = pl,
        well_pos = wp,
        time     = time_pts,
        od       = od,
        Chemical = "CONTROL",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }

    row_of <- function(wp) ceiling(as.integer(wp) / 12L)

    for (i in seq_along(cmpd_pos)) {
      wp     <- cmpd_pos[i]
      cmpd   <- cmpd_names[i]
      is_hit <- cmpd %in% hits_pl
      k_val  <- if (is_hit) runif(1, 0.05, 0.35) else runif(1, 0.65, 0.95)
      od     <- logistic_od(time_pts, k = k_val, t_mid = 18 + rnorm(1, 0, 2)) +
                row_of(wp) * 0.003
      all_rows[[idx]] <- data.frame(
        run      = run_id,
        plate    = pl,
        well_pos = wp,
        time     = time_pts,
        od       = od,
        Chemical = cmpd,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }

    example_screen[[run_id]] <- do.call(rbind, all_rows)
  }
}

dir.create("data", showWarnings = FALSE)
save(example_screen, file = "data/example_screen.rda", compress = "bzip2")

cat("Saved data/example_screen.rda\n")
cat(sprintf("Runs: %d  |  Rows/run: %d  |  Wells/run: %d\n",
            length(example_screen),
            nrow(example_screen[[1]]),
            length(unique(example_screen[[1]]$well_pos))))
