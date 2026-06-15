## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width = 7,
  fig.height = 4
)


## ----quickstart, eval = FALSE-------------------------------------------------
# library(gcscreen)
# data(example_screen)
# 
# res <- run_screen(
#   platerun_list = example_screen,
#   time_max      = 60,
#   fdr_threshold = 0.05,
#   auc_rel_max   = 0.99
# )
# 
# # Confirmed hits (ANOVA ∩ LMM, FDR < 0.05)
# res$hits[, c("Chemical", "AUC_rel", "FDR", "FDR_lmm")]


## ----layout-------------------------------------------------------------------
library(gcscreen)

# Map well positions to row/column indices
well_to_rowcol(c(1L, 12L, 13L, 96L))


## ----auc----------------------------------------------------------------------
data(example_screen)
run_df <- example_screen[["placa_2_rep_A"]]

# Select one well, sort by time, baseline-correct, compute AUC
w1 <- run_df[run_df$well_pos == 13L, ]   # a control well (col 1, row 2)
w1 <- w1[order(w1$time), ]
od_corr <- w1$od - w1$od[1L]
od_corr[od_corr < 0] <- 0

auc_val <- auc_trapz(w1$time, od_corr)
cat("AUC for well 13:", round(auc_val, 3), "\n")


## ----bscore-------------------------------------------------------------------
# Simulate a plate-run AUC vector with a row gradient
set.seed(42)
wpos    <- 1:96
rc      <- well_to_rowcol(wpos)
auc_raw <- rnorm(96, mean = 40, sd = 2) + rc$row_idx * 1.5

auc_corr <- bscore_plate(wpos, auc_raw)

# Row-AUC correlation in compound wells before and after correction
cmpd_idx <- which(rc$col_idx %in% 2:11)
cat("Pearson r(row, AUC) before B-score:",
    round(cor(rc$row_idx[cmpd_idx], auc_raw[cmpd_idx]),  3), "\n")
cat("Pearson r(row, AUC) after  B-score:",
    round(cor(rc$row_idx[cmpd_idx], auc_corr[cmpd_idx]), 3), "\n")


## ----qc-----------------------------------------------------------------------
# Compute AUC for all control wells of one plate-run
ctrl_df <- run_df[run_df$Chemical == "CONTROL", ]
ctrl_auc <- tapply(seq_len(nrow(ctrl_df)), ctrl_df$well_pos, function(idx) {
  sub <- ctrl_df[idx, ]
  sub <- sub[order(sub$time), ]
  od  <- sub$od - sub$od[1L]; od[od < 0] <- 0
  auc_trapz(sub$time, od)
})

qc <- plate_qc(as.numeric(ctrl_auc), cv_warn = 20)
cat(sprintf("CV: %.1f%%  |  flag: %s\n", qc$cv, qc$cv_flag))


## ----normalize----------------------------------------------------------------
# Normalise compound AUC to the control median
ctrl_vec <- as.numeric(ctrl_auc)
cmpd_vec <- c(2.1, 18.4, 35.0, 38.9)   # example values
auc_rel  <- normalize_controls(cmpd_vec, ctrl_vec)
cat("AUC_rel:", round(auc_rel, 3), "\n")


## ----stats-single-------------------------------------------------------------
# ANOVA: 3 compound replicates vs 48 control replicates (3 plate-runs × 16)
ctrl_all <- rnorm(48, mean = 40, sd = 3)
cmpd_obs <- c(4.2, 3.8, 4.5)            # a strong inhibitor
pval <- screen_anova(cmpd_obs, ctrl_all)
cat("ANOVA p-value:", format(pval, digits = 3), "\n")


## ----full-screen, eval = FALSE------------------------------------------------
# data(example_screen)
# res <- run_screen(example_screen)
# 
# # Inspect QC table
# head(res$qc)
# 
# # Top hits by ANOVA FDR
# head(res$anova_results)
# 
# # Confirmed hits (ANOVA ∩ LMM)
# res$hits


## ----hit-calling, eval = FALSE------------------------------------------------
# hits <- call_hits(
#   anova_res     = res$anova_results,
#   lmm_res       = res$lmm_results,
#   fdr_threshold = 0.05,
#   auc_rel_max   = 0.99,
#   exclude       = "Chicago sky blue 6B"   # optical artefact
# )

