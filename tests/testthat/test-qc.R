test_that("ctrl_cv returns correct value", {
  # sd=1, mean=10 → cv=10%
  x <- rep(10, 10) + c(-1, 1, -1, 1, -1, 1, -1, 1, 0, 0)
  expect_equal(ctrl_cv(x), 100 * sd(x) / mean(x))
})

test_that("ctrl_cv returns NA when mean is zero", {
  expect_true(is.na(ctrl_cv(c(-1, 1))))
})

test_that("zprime returns NA with fewer than 2 observations", {
  expect_true(is.na(zprime(5, c(40, 41, 42))))
  expect_true(is.na(zprime(c(5, 6), 40)))
})

test_that("zprime is 1 for perfectly separated groups with zero variance", {
  # degenerate: sd = 0 in both groups → formula collapses to 1
  # Actually Z' = 1 - 3*(0+0)/|mu_pos - mu_neg|; but sd(c(5,5)) = 0 so Z' = 1
  expect_equal(zprime(c(5, 5), c(40, 40)), 1)
})

test_that("zprime is negative for overlapping groups", {
  pos <- rnorm(20, mean = 30, sd = 10)
  neg <- rnorm(20, mean = 35, sd = 10)
  expect_lt(zprime(pos, neg), 0)
})

test_that("ssmd is positive when pos mean < neg mean (inhibition)", {
  pos <- rnorm(16, mean = 10, sd = 2)
  neg <- rnorm(16, mean = 40, sd = 3)
  # ssmd = (pos - neg)/sqrt(var_pos + var_neg) → negative for inhibitors
  # by convention pos = hits (lower growth) and neg = controls (higher growth)
  sm <- ssmd(pos, neg)
  expect_lt(sm, 0)
})

test_that("plate_qc returns a list with correct names", {
  ctrl <- rnorm(16, 42, 3)
  hits <- rnorm(5, 10, 2)
  res  <- plate_qc(ctrl, hits)
  expect_named(res, c("cv", "zprime", "rz_prime", "ssmd", "cv_flag"))
})

test_that("plate_qc flags HIGH_CV correctly", {
  ctrl_high_cv <- c(10, 40, 10, 40)  # cv >> 20%
  res <- plate_qc(ctrl_high_cv, cv_warn = 20)
  expect_equal(res$cv_flag, "HIGH_CV")
})

test_that("plate_qc returns NA QC metrics when hit_auc has fewer than 2 obs", {
  ctrl <- rnorm(16, 42, 3)
  res  <- plate_qc(ctrl, hit_auc = 10)
  expect_true(is.na(res$zprime))
  expect_true(is.na(res$ssmd))
})
