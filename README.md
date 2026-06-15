# gcscreen — Kinetic Growth-Curve Analysis for 96-Well Phenotypic Screens

[![R-CMD-check](https://github.com/mparmol/gcscreen/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mparmol/gcscreen/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20702376.svg)](https://doi.org/10.5281/zenodo.20702376)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

`gcscreen` is an R package implementing a complete, reproducible
pipeline for phenotypic screening using OD600 kinetic growth curves
in 96-well plates. It provides trapezoidal AUC calculation, B-score
positional-bias correction (median polish), control normalisation,
plate-quality metrics (CV, Z'-factor, SSMD, robust Z'-factor),
per-compound ANOVA and linear mixed-model testing (with plate as
random effect), and Benjamini–Hochberg FDR-based hit calling.
Designed for *Saccharomyces cerevisiae* high-throughput screens
but applicable to any kinetic OD-based 96-well assay.

## Why this package

Plate-based phenotypic screens with kinetic OD600 readout are
common in chemical biology, but their validity hinges on the
control of plate-level technical variability — a problem with
standard solutions in single-endpoint high-throughput screening
that are unevenly adopted in growth-curve work. `gcscreen`
implements the recommended best practices (B-score correction,
median-of-controls normalisation, FDR-controlled per-plate ANOVA
and library-wide linear mixed model with plate as random effect)
in a single coherent pipeline.

The companion paper benchmarks the pipeline against five
alternative normalisation strategies, validates its operating
characteristics by Monte Carlo simulation, and confirms — through
functional principal component analysis and parametric logistic
fitting — that the area under the growth curve is the appropriate
summary statistic for this assay format.

## Installation

```r
# Install from GitHub
remotes::install_github("mparmol/gcscreen")

# Or, after cloning locally:
devtools::install("path/to/gcscreen")
```

R >= 4.1.0 required. Imports `lmerTest` (>= 3.1.0) and base R
`stats`. Suggests `testthat`, `knitr`, `rmarkdown`, `ggplot2` for
testing, vignette building and visualisation.

## Quick start

```r
library(gcscreen)

# Built-in synthetic example: 9 plate-runs, 12 planted hits, 228 non-hits
data(example_screen)
str(example_screen)

# Run the complete pipeline in one line
result <- run_screen(example_screen)

# Inspect results
result$hits        # core hit list (ANOVA ∩ LMM)
result$qc          # plate quality metrics
result$auc         # per-well AUC after B-score + normalisation
```

For a full walk-through, see the package vignette:

```r
vignette("gcscreen-workflow")
```

## Exported functions

| Function                      | Purpose                                                  |
| ----------------------------- | -------------------------------------------------------- |
| `auc_trapz(time, od)`         | Trapezoidal AUC of an OD600 trace                        |
| `baseline_correct(time, od)`  | Subtract the t = 0 baseline                              |
| `bscore_plate(well_pos, auc)` | B-score (median polish) on the 8 × 10 compound submatrix |
| `normalize_controls(...)`     | Divide compound AUCs by the median of in-plate controls  |
| `plate_qc(ctrl, hit)`         | All four plate-quality metrics in one call               |
| `zprime`, `rz_prime`, `ssmd`  | Individual quality metrics                               |
| `ctrl_cv(ctrl)`               | Coefficient of variation of negative controls            |
| `screen_anova(...)`           | One-way ANOVA per compound + BH-FDR                      |
| `screen_lmm(...)`             | Linear mixed model per compound + BH-FDR                 |
| `run_anova_screen(df)`        | Apply ANOVA across the whole library                     |
| `run_lmm_screen(df)`          | Apply LMM across the whole library                       |
| `call_hits(anova_res, lmm_res)` | Convergence-based hit calling (ANOVA ∩ LMM)            |
| `run_screen(df)`              | Single-call wrapper for the complete pipeline            |

See `?function_name` for full documentation. All functions are
documented with roxygen2 and have unit-test coverage.

## Citation

If you use `gcscreen`, please cite both the package and the
companion paper:

- **Package**: Parras-Moltó M (2026). *gcscreen: Kinetic growth-curve
  analysis for 96-well phenotypic screens*. R package version
  0.1.0. doi:10.5281/zenodo.20702376.
- **Paper**: Parras-Moltó M, García-Ríos E (2026). *Triangulating
  compound effects in 96-well growth-curve phenotypic screens: a
  benchmarked, simulation-validated pipeline*. PLOS Computational
  Biology (in submission). doi:10.5281/zenodo.20702473.

A `CITATION.cff` file is provided for automatic citation tools
(GitHub "Cite this repository" button, Zotero, etc.).

## Reproducibility

This package was developed and tested with R 4.5.x. The complete
session info under which the companion paper analyses were run is
deposited in the [paper repository][paper-repo] under
`environment/sessionInfo.txt`.

[paper-repo]: https://github.com/mparmol/prestwick-yeast-screen

## Contributing

Issues and pull requests welcome. Please run
`devtools::check()` and ensure `R CMD check --as-cran` passes
(0 errors / 0 warnings) before opening a PR.

## Authors

- **Marcos Parras-Moltó** (creator and maintainer) —
  mparmol@gmail.com — ORCID
  [0000-0003-0529-627X](https://orcid.org/0000-0003-0529-627X) —
  Leitat Technological Center, Applied Microbiology and
  Biotechnologies, Barcelona, Spain.
- **Estéfani García-Ríos** (co-author of the companion paper) —
  egarcia@icvv.es — ORCID
  [0000-0001-9028-055X](https://orcid.org/0000-0001-9028-055X) —
  Instituto de Ciencias de la Vid y del Vino (CSIC, Universidad
  de La Rioja, Gobierno de La Rioja), Logroño 26007, Spain.

## License

MIT © 2026 Marcos Parras-Moltó. See `LICENSE` for full text.
