
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DeSurv

<!-- badges: start -->

<!-- badges: end -->

DeSurv couples nonnegative matrix factorisation (NMF) with a Cox
proportional hazards layer so that latent expression programs are
learned with survival in mind. The package currently provides:

- Fast C++ implementations of the alternating updates used to fit the
  latent factors, including elastic-net penalised Cox coefficients.
- Automatic multi-start initialisation and warm-starts across a grid of
  supervision strengths.
- Cross-validation utilities that summarise concordance indices across
  folds/inits and refit the winning model on the full cohort.
- Prediction helpers for downstream survival pipelines and lightweight
  feature pre-processing (e.g. gene filtering).

## Installation

The package is not on CRAN yet. Install the development version from
GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ayoung31/DeSurv")
```

## Getting started

``` r
library(DeSurv)

set.seed(123)
p <- 50; n <- 120; k <- 3
X <- matrix(rexp(p * n), nrow = p, ncol = n)
rownames(X) <- paste0("gene", seq_len(p))
y <- rexp(n) + 0.1
d <- rbinom(n, 1, 0.7)

dd <- desurv_data(X, y, d, k = k)

fit <- desurv_fit(
  dd,
  alpha   = 0.5,
  lambda  = 1,
  nu      = 0.5,
  lambdaW = 0.1,
  lambdaH = 0.1,
  ninit   = 5,
  imaxit  = 15,
  maxit   = 200,
  verbose = FALSE
)

head(predict(fit, type = "risk"))
```

## Cross-validation

``` r
cv_fit <- desurv_cv(
  dd,
  k_grid       = 2:3,
  alpha_grid   = c(0.25, 0.5, 0.75),
  lambda_grid  = c(0.5, 1),
  nu_grid      = c(0.25, 0.5),
  lambdaW_grid = c(0.05, 0.1),
  lambdaH_grid = c(0.05, 0.1),
  n_starts     = 2,
  nfolds       = 3,
  tol          = 1e-5,
  maxit        = 300,
  verbose      = FALSE
)

cv_fit$best
cv_fit$fit$cindex
```

See the documentation for `desurv_fit()`, `desurv_alpha_warmstart()`,
and `desurv_cv()` for additional arguments such as custom initial
values, parallelisation, or stricter convergence thresholds.
