test_that("screen_anova returns low p-value for clearly different groups", {
  cmpd <- rnorm(3, mean = 5,  sd = 1)
  ctrl <- rnorm(48, mean = 40, sd = 3)
  pval <- screen_anova(cmpd, ctrl)
  expect_lt(pval, 0.001)
})

test_that("screen_anova returns high p-value for similar groups", {
  set.seed(99)
  cmpd <- rnorm(3, mean = 40, sd = 3)
  ctrl <- rnorm(48, mean = 40, sd = 3)
  pval <- screen_anova(cmpd, ctrl)
  expect_gt(pval, 0.05)
})

test_that("screen_anova returns NA with too few observations", {
  expect_true(is.na(screen_anova(5, c(40, 41))))   # 1 compound obs
  expect_true(is.na(screen_anova(c(5, 6), 40)))    # 1 control obs
})

test_that("screen_lmm uses ANOVA fallback for single-plate compound", {
  set.seed(7)
  df <- data.frame(
    auc      = c(rnorm(3, 15, 2), rnorm(16, 42, 3)),
    Chemical = c(rep("DrugX", 3), rep("CONTROL", 16)),
    plate    = 2L
  )
  res <- screen_lmm(df, compound_name = "DrugX")
  expect_true(res$method %in% c("ANOVA", "ANOVA_fallback"))
  expect_lt(res$p_value, 0.05)
})

test_that("screen_lmm uses LMM for multi-plate compound", {
  set.seed(8)
  df <- data.frame(
    auc      = c(rnorm(3, 15, 2), rnorm(48, 42, 3)),
    Chemical = c(rep("DrugY", 3), rep("CONTROL", 48)),
    plate    = c(2L, 3L, 5L, sample(c(2L, 3L, 5L), 48, replace = TRUE))
  )
  res <- screen_lmm(df, compound_name = "DrugY")
  expect_equal(res$method, "LMM")
  expect_lt(res$p_value, 0.05)
})

test_that("call_hits returns intersection of ANOVA and LMM", {
  anova_res <- data.frame(
    Chemical = c("A", "B", "C", "D"),
    Plate    = 2L, n_reps = 3L,
    AUC_rel  = c(0.1, 0.5, 0.2, 0.8),
    p_anova  = c(0.001, 0.002, 0.003, 0.004),
    FDR      = c(0.004, 0.008, 0.012, 0.016),
    stringsAsFactors = FALSE
  )
  lmm_res <- data.frame(
    Chemical = c("A", "B", "C", "D"),
    FDR      = c(0.003, 0.06, 0.01, 0.02),   # B fails LMM FDR
    stringsAsFactors = FALSE
  )
  hits <- call_hits(anova_res, lmm_res, fdr_threshold = 0.05)
  expect_false("B" %in% hits$Chemical)   # B fails LMM
  expect_true("A" %in% hits$Chemical)
  expect_true("C" %in% hits$Chemical)
})

test_that("call_hits respects exclude argument", {
  anova_res <- data.frame(
    Chemical = "DrugA", Plate = 2L, n_reps = 3L,
    AUC_rel = 0.1, p_anova = 0.001, FDR = 0.002,
    stringsAsFactors = FALSE
  )
  lmm_res <- data.frame(Chemical = "DrugA", FDR = 0.002,
                         stringsAsFactors = FALSE)
  hits <- call_hits(anova_res, lmm_res, exclude = "DrugA")
  expect_equal(nrow(hits), 0L)
})

test_that("call_hits filters by auc_rel_max", {
  anova_res <- data.frame(
    Chemical = c("Inhibitor", "Activator"),
    Plate = 2L, n_reps = 3L,
    AUC_rel = c(0.2, 1.3),
    p_anova = c(0.001, 0.001),
    FDR     = c(0.002, 0.002),
    stringsAsFactors = FALSE
  )
  lmm_res <- data.frame(
    Chemical = c("Inhibitor", "Activator"),
    FDR      = c(0.002, 0.002),
    stringsAsFactors = FALSE
  )
  hits <- call_hits(anova_res, lmm_res)
  expect_true("Inhibitor" %in% hits$Chemical)
  expect_false("Activator" %in% hits$Chemical)
})
