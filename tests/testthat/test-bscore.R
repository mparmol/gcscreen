test_that("bscore_plate returns same length as input", {
  set.seed(1)
  auc <- rnorm(96, 40, 4)
  out <- bscore_plate(1:96, auc)
  expect_equal(length(out), 96L)
})

test_that("bscore_plate does not modify control wells", {
  set.seed(2)
  auc <- rnorm(96, 40, 4)
  out <- bscore_plate(1:96, auc)
  rc  <- well_to_rowcol(1:96)
  ctrl_idx <- which(rc$col_idx %in% c(1L, 12L))
  expect_equal(out[ctrl_idx], auc[ctrl_idx])
})

test_that("bscore_plate reduces row-gradient in compound wells", {
  set.seed(3)
  wpos <- 1:96
  rc   <- well_to_rowcol(wpos)
  # Pure row gradient: each row adds 2 AUC units
  auc  <- 40 + rc$row_idx * 2 + rnorm(96, 0, 0.1)
  out  <- bscore_plate(wpos, auc)

  cmpd_idx <- which(rc$col_idx %in% 2:11)
  # After correction the row gradient should be near zero
  cor_before <- cor(rc$row_idx[cmpd_idx], auc[cmpd_idx])
  cor_after  <- cor(rc$row_idx[cmpd_idx], out[cmpd_idx])
  expect_gt(abs(cor_before), 0.9)
  expect_lt(abs(cor_after),  0.2)
})

test_that("bscore_plate warns and returns input when too few compound wells", {
  # Use only control-column positions (cols 1 and 12): 0 compound wells
  ctrl_wells <- c(1L, 13L, 25L, 37L, 49L, 61L, 73L, 85L,
                  12L, 24L, 36L, 48L, 60L, 72L, 84L, 96L)
  auc <- rnorm(16L, 40, 4)
  expect_warning(out <- bscore_plate(ctrl_wells, auc))
  expect_equal(out, auc)
})
