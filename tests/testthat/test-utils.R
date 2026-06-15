test_that("well_to_rowcol maps corners correctly", {
  expect_equal(well_to_rowcol(1L),  data.frame(row_idx = 1L, col_idx = 1L))
  expect_equal(well_to_rowcol(12L), data.frame(row_idx = 1L, col_idx = 12L))
  expect_equal(well_to_rowcol(13L), data.frame(row_idx = 2L, col_idx = 1L))
  expect_equal(well_to_rowcol(96L), data.frame(row_idx = 8L, col_idx = 12L))
})

test_that("well_to_rowcol rejects out-of-range values", {
  expect_error(well_to_rowcol(0L))
  expect_error(well_to_rowcol(97L))
})

test_that("auc_trapz returns correct area for linear ramp", {
  # triangle: base = 10, height = 1 → area = 5
  t  <- c(0, 10)
  od <- c(0, 1)
  expect_equal(auc_trapz(t, od), 5)
})

test_that("auc_trapz returns NA for single point", {
  expect_true(is.na(auc_trapz(0, 0.1)))
})

test_that("auc_trapz errors on length mismatch", {
  expect_error(auc_trapz(1:3, 1:2))
})

test_that("baseline_correct subtracts first row and floors at zero", {
  m <- matrix(c(0.1, 0.2, 0.3, 0.05, 0.15, 0.25), nrow = 3, ncol = 2)
  out <- baseline_correct(m)
  expect_equal(out[1L, ], c(0, 0))
  expect_equal(out[2L, ], c(0.1, 0.1))
  expect_true(all(out >= 0))
})
