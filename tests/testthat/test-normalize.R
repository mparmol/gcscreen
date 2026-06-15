test_that("normalize_controls divides by median control", {
  ctrl <- c(40, 42, 38, 44)   # median = 41
  cmpd <- c(20, 41, 82)
  out  <- normalize_controls(cmpd, ctrl)
  med  <- median(ctrl)
  expect_equal(out, cmpd / med)
})

test_that("normalize_controls returns NA when control median is zero", {
  expect_warning(out <- normalize_controls(c(1, 2), c(0, 0)))
  expect_true(all(is.na(out)))
})

test_that("normalize_controls handles NA in controls", {
  ctrl <- c(40, NA, 42, 38)
  cmpd <- c(20, 40)
  out  <- normalize_controls(cmpd, ctrl)
  expect_equal(out, cmpd / median(ctrl, na.rm = TRUE))
})
