# gcscreen 0.1.0 (2026-05-XX)

Initial public release accompanying the manuscript
*Triangulating compound effects in 96-well growth-curve phenotypic
screens: a benchmarked, simulation-validated pipeline combining
within-plate B-score correction, median-of-controls normalisation
and linear mixed modelling* (Parras-Moltó & García-Ríos, 2026).

## New features

- Trapezoidal AUC computation (`auc_trapz()`) and baseline
  correction (`baseline_correct()`) for OD600 time series.
- Within-plate B-score positional-bias correction by iterative
  median polish on the 8 × 10 compound submatrix
  (`bscore_plate()`).
- Median-of-controls normalisation (`normalize_controls()`)
  replacing single-reference-well alternatives.
- Plate-level quality metrics: CV, Z'-factor, SSMD and robust
  Z'-factor (`plate_qc()`, `zprime()`, `rz_prime()`, `ssmd()`,
  `ctrl_cv()`).
- Per-compound statistical testing: one-way ANOVA against in-plate
  controls (`screen_anova()`) and linear mixed model with plate
  as random effect (`screen_lmm()`, via `lmerTest`).
- Library-wide BH-FDR correction and convergence-based hit calling
  ANOVA ∩ LMM (`call_hits()`).
- Single-call complete pipeline wrapper (`run_screen()`).
- Built-in synthetic example dataset (`data(example_screen)`):
  9 plate-runs, 12 planted hits at known effect sizes, 228 non-hit
  compounds.
- Comprehensive test suite (47 unit tests, `testthat` edition 3)
  covering all exported functions including edge cases.
- Vignette `gcscreen-workflow` demonstrating the full analysis
  pipeline from raw OD600 trajectories to hit list.
- `R CMD check --as-cran`: 0 errors / 0 warnings / 1 note (future
  timestamps inherent to release date).

## Documentation

- Roxygen2-generated `man/*.Rd` for all 16 exported functions.
- README with installation, quick-start example, function
  reference table and citation block.

## License

MIT © 2026 Marcos Parras-Moltó.
