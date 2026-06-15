# Shared helpers for gcscreen tests

# Minimal valid plate-run: 96 wells, 16 controls (cols 1 & 12), 80 compounds
make_platerun <- function(plate_id = 2L, run_id = "placa_2_rep_A",
                          n_hits = 0L, seed = 1L) {
  set.seed(seed)
  time_pts <- seq(0, 60, by = 4)

  ctrl_pos <- c(1:8, 89:96)
  cmpd_pos <- 9:88
  n_cmpd   <- length(cmpd_pos)

  rows <- NULL

  for (wp in ctrl_pos) {
    od <- pmax(plogis((time_pts - 18) / 5) * 0.8 + rnorm(length(time_pts), 0, 0.005), 0)
    rows <- rbind(rows, data.frame(run = run_id, plate = plate_id,
                                   well_pos = wp, time = time_pts,
                                   od = od, Chemical = "CONTROL",
                                   stringsAsFactors = FALSE))
  }

  for (i in seq_len(n_cmpd)) {
    wp      <- cmpd_pos[i]
    is_hit  <- i <= n_hits
    k_val   <- if (is_hit) 0.15 else 0.80
    chem    <- if (is_hit) paste0("Hit_", i) else paste0("Cmpd_", i)
    od <- pmax(plogis((time_pts - 18) / 5) * k_val + rnorm(length(time_pts), 0, 0.005), 0)
    rows <- rbind(rows, data.frame(run = run_id, plate = plate_id,
                                   well_pos = wp, time = time_pts,
                                   od = od, Chemical = chem,
                                   stringsAsFactors = FALSE))
  }
  rows
}
