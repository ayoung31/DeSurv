# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Development Commands

```r
# Load package for interactive development (no reinstall needed)
devtools::load_all()

# Regenerate documentation after editing roxygen comments
devtools::document()

# Run all tests
devtools::test()

# Run a single test file
testthat::test_file("tests/testthat/test-fit-workflow.R")

# Full package check (build, test, lint)
devtools::check()

# Recompile C++ code quickly during development
Rcpp::sourceCpp('src/functions.cpp')

# Regenerate Rcpp bindings (never edit RcppExports.* manually)
Rcpp::compileAttributes()

# Check test coverage
covr::report()
```

## Architecture Overview

DeSurv implements survival-driven nonnegative matrix factorization (NMF), coupling latent factor learning with a Cox proportional hazards model.

**Core Model**: Decomposes nonnegative matrix X (p×n, genes×samples) into W (p×k) and H (k×n), fitting a Cox model with linear predictor η = (X'W)β where β ∈ Rᵏ.

### Optimization Flow

The model uses alternating block updates in `optimize_loss_cpp()`:

1. **H update** (`update_H_cpp`): Multiplicative NMF update for the coefficient matrix
2. **W update** (`update_W_damped_backtrack`): Damped multiplicative update with backtracking line search; combines NMF gradient with Cox gradient (scaled to balance magnitudes)
3. **β update** (`update_beta_cpp`): Coordinate descent with elastic-net penalty on standardized Z = X'W; uses log-sum-exp suffix statistics for numerical stability

Convergence is based on cosine distance between successive W and H iterates.

### Key Components

| File | Purpose |
|------|---------|
| `R/desurv_fit.R` | Main fitting function; multi-start initialization → full optimization |
| `R/init.R` | Short warm-up runs to find good starting (W, H, β) |
| `R/desurv_cv.R` | K-fold CV with warmstart or cold engine; model selection via max or 1-SE rule |
| `R/cv_warm.R`, `R/cv_cold.R` | Internal CV engines (warmstart reuses α paths; cold refits independently) |
| `R/desurv_alpha_warmstart.R` | Warm-started optimization across supervision weights (α grid) |
| `R/desurv_bayesopt.R`, `R/desurv_cv_bayesopt.R` | Bayesian optimization for hyperparameter tuning (requires DiceKriging, lhs) |
| `R/predict_methods.R` | S3 `predict()` and `coef()` methods |
| `R/preprocess_data.R` | Gene filtering, rank/quantile normalization |
| `R/optimizer_wrappers.R` | R wrapper `.run_optimize_loss()` that calls C++ |
| `src/functions.cpp` | All performance-critical routines (Rcpp/RcppArmadillo) |

### Data Flow

1. `desurv_data()` validates and packages X, y, d, k into a data object
2. `desurv_fit()` runs multi-start initialization (`init()`) then full optimization
3. Prediction via `predict.desurv_fit()` returns risk scores or linear predictors

### C++ Implementation Details

`src/functions.cpp` contains:

- **`cox_suffix_stats_full`**: Computes gradient g_z, diagonal Hessian approx w, and partial log-likelihood using log-sum-exp with suffix max stabilization. Handles tied event times via `gstart` index.
- **`cdfit_cox_dh_one_lambda`**: Coordinate descent for elastic-net Cox regression with proximal (soft-threshold) updates
- **`calc_loss_cpp`**: Combined objective: `(1-α)(NMF_loss + λH·pen_H) - α·surv_loss + λW·pen_W`

### Hyperparameters

| Parameter | Range | Description |
|-----------|-------|-------------|
| **alpha** | [0,1] | Supervision strength (0 = pure NMF, 1 = pure Cox) |
| **nu** | [0,1] | Elastic-net mixing for β (0 = ridge, 1 = lasso) |
| **lambda** | ≥0 | Global penalty scale for β |
| **lambdaW** | ≥0 | L2 penalty on W |
| **lambdaH** | ≥0 | L2 penalty on H |

### Cross-Validation Engines

- **warmstart**: Runs `desurv_alpha_warmstart()` per fold; reuses (W, H, β) along α path
- **cold**: Calls `desurv_fit()` independently for each α; `n_starts` maps to `ninit`

Model selection uses either `"max"` (best mean C-index) or `"1se"` (simplest within 1 SE of best) with tie-breaking: smallest k → smallest α → largest λ → smallest ν → smallest λW → smallest λH.

## Coding Conventions

- **R**: tidyverse style, 2-space indent, snake_case, S3 methods as `generic.class`
- **C++**: camelCase functions, `const auto&` where possible, guard exports with `// [[Rcpp::export]]` docs
- Input validation via helpers in `R/validate_input.R`
- Run `styler::style_pkg()` before large refactors

## Testing

Tests in `tests/testthat/test-*.R`. Key test files:
- `test-fit-workflow.R`: End-to-end fitting
- `test-cv-folds.R`: Cross-validation fold logic
- `test-bayesopt.R`, `test-cv-bayesopt.R`: Bayesian optimization

Target coverage ≥85% for statistical cores.

## Dependencies

**Required**: Rcpp, RcppArmadillo, RcppEigen, survival, cvwrapr, preprocessCore

**Optional** (for Bayesian optimization): DiceKriging, lhs
