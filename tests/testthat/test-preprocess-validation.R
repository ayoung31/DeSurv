# Tests for preprocessing validation fixes
# Covers Findings 1-6 from the preprocessing & prediction validation review

# Helper to create valid test data
make_valid_preprocess_data <- function(n = 20, p = 10, seed = 42) {
  set.seed(seed)
  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  rownames(X) <- paste0("gene", seq_len(p))
  colnames(X) <- paste0("sample", seq_len(n))
  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.6)
  if (sum(d) < 3) d[1:3] <- 1
  dataset <- rep(c("A", "B"), each = n / 2)
  list(X = X, y = y, d = d, dataset = dataset)
}

# =============================================================================
# Finding 1: preprocess_data() validates y, d, and dataset
# =============================================================================

test_that("preprocess_data() errors on non-numeric y", {
  data <- make_valid_preprocess_data()
  data$y <- as.character(data$y)
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`y` must be a numeric vector"
  )
})

test_that("preprocess_data() errors on NA in y", {
  data <- make_valid_preprocess_data()
  data$y[5] <- NA
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`y` cannot contain NA"
  )
})

test_that("preprocess_data() errors on Inf in y", {
  data <- make_valid_preprocess_data()

  # Test Inf
  data$y[5] <- Inf
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`y` cannot contain Inf or NaN"
  )
})

test_that("preprocess_data() errors on NaN in y", {
  data <- make_valid_preprocess_data()

  # Test NaN - in R, is.na(NaN) == TRUE, so it gets caught by NA check
  data$y[5] <- NaN
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`y` cannot contain NA"  # NaN is considered NA in R
  )
})

test_that("preprocess_data() errors on non-positive y", {
  data <- make_valid_preprocess_data()
  data$y[5] <- 0
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`y` \\(survival times\\) must be positive"
  )

  data <- make_valid_preprocess_data()
  data$y[5] <- -5
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`y` \\(survival times\\) must be positive"
  )
})

test_that("preprocess_data() errors on NA in d", {
  data <- make_valid_preprocess_data()
  data$d[5] <- NA
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`d` cannot contain NA"
  )
})

test_that("preprocess_data() errors on invalid d values", {
  data <- make_valid_preprocess_data()
  data$d[5] <- 2
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`d` must contain only 0 \\(censored\\) and 1 \\(event\\)"
  )

  data <- make_valid_preprocess_data()
  data$d[5] <- -1
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`d` must contain only 0 \\(censored\\) and 1 \\(event\\)"
  )

  data <- make_valid_preprocess_data()
  data$d[5] <- 0.5
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`d` must contain only 0 \\(censored\\) and 1 \\(event\\)"
  )
})

test_that("preprocess_data() errors on NA in dataset", {
  data <- make_valid_preprocess_data()
  data$dataset[5] <- NA
  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    ),
    "`dataset` cannot contain NA"
  )
})

# =============================================================================
# Finding 2: preprocess_data() handles Inf/NaN in X, not just NA
# =============================================================================

test_that("preprocess_data() replaces Inf values in X with 0", {
  data <- make_valid_preprocess_data()
  data$X[1, 1] <- Inf
  data$X[2, 2] <- -Inf

  result <- preprocess_data(
    X = data$X, y = data$y, d = data$d, dataset = data$dataset,
    method_trans_train = "none", verbose = FALSE
  )

  # Result should not contain Inf
  expect_true(all(is.finite(result$ex)))
})

test_that("preprocess_data() replaces NaN values in X with 0", {
  data <- make_valid_preprocess_data()
  data$X[1, 1] <- NaN

  result <- preprocess_data(
    X = data$X, y = data$y, d = data$d, dataset = data$dataset,
    method_trans_train = "none", verbose = FALSE
  )

  # Result should not contain NaN
  expect_true(all(!is.nan(result$ex)))
})

test_that("preprocess_data() replaces NA values in X with 0", {
  data <- make_valid_preprocess_data()
  data$X[1, 1] <- NA

  result <- preprocess_data(
    X = data$X, y = data$y, d = data$d, dataset = data$dataset,
    method_trans_train = "none", verbose = FALSE
  )

  # Result should not contain NA
  expect_true(all(!is.na(result$ex)))
})

test_that("preprocess_data() emits message about non-finite values when verbose", {
  data <- make_valid_preprocess_data()
  data$X[1, 1] <- Inf
  data$X[2, 2] <- NA

  expect_message(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = TRUE
    ),
    "Replaced.*non-finite values"
  )
})

# =============================================================================
# Finding 6: Explicit genes list warns when genes are dropped
# =============================================================================

test_that("preprocess_data() warns when explicit genes are dropped", {
  data <- make_valid_preprocess_data()
  # Request genes, some of which don't exist
  genes <- c("gene1", "gene2", "NONEXISTENT_GENE", "ANOTHER_FAKE_GENE")

  expect_warning(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      genes = genes,
      method_trans_train = "none", verbose = FALSE
    ),
    "genes not found in X and were dropped"
  )
})

test_that("preprocess_data() does not warn when all explicit genes are present", {
  data <- make_valid_preprocess_data()
  genes <- c("gene1", "gene2", "gene3")

  expect_no_warning(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      genes = genes,
      method_trans_train = "none", verbose = FALSE
    )
  )
})

test_that("preprocess_data() errors when NO explicit genes are found", {
  data <- make_valid_preprocess_data()
  genes <- c("FAKE1", "FAKE2", "FAKE3")

  expect_error(
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      genes = genes,
      method_trans_train = "none", verbose = FALSE
    ),
    "None of the supplied `genes` are present"
  )
})

# =============================================================================
# Finding 3: preprocess_info timing - predict works with quant normalization
# =============================================================================

test_that("desurv_cv with preprocess=TRUE and quant method creates usable preprocess_info", {
  skip_if_not_installed("preprocessCore")

  # Create small dataset for CV
  data <- make_valid_preprocess_data(n = 16, p = 20)

  # Run CV with preprocessing and quantile normalization
  cv_result <- suppressWarnings(desurv_cv(
    X = data$X,
    y = data$y,
    d = data$d,
    dataset = data$dataset,
    k_grid = 2,
    alpha_grid = 0.5,
    lambda_grid = 0.1,
    nu_grid = 0.5,
    lambdaW_grid = 0,
    lambdaH_grid = 0,
    preprocess = TRUE,
    method_trans_train = "quant",
    nfolds = 2,
    tol = 1e-3,
    maxit = 30,
    verbose = FALSE
  ))

  # The fit should have complete preprocess metadata
  fit <- cv_result$fit
  expect_true(!is.null(fit$preprocess))
  expect_equal(fit$preprocess$method_trans_train, "quant")
  expect_true(!is.null(fit$preprocess$transform_target))
  expect_true(!is.null(fit$preprocess$transform_target$values))

  # Prediction should work without "Quantile normalization target is missing" error
  newX <- data$X[, 1:5, drop = FALSE]
  expect_no_error({
    preds <- predict(fit, newdata = newX)
  })
  expect_length(preds, 5)
})

# =============================================================================
# Finding 4: desurv_data warning on preprocessing bypass
# =============================================================================

test_that("predict warns when desurv_data bypasses preprocessing transforms", {
  data <- make_valid_preprocess_data()

  # Create a fit with preprocessing metadata
  fit <- suppressWarnings(desurv_fit(
    X = data$X, y = data$y, d = data$d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 30, ninit = 1,
    verbose = FALSE,
    preprocess_info = list(
      method_trans_train = "rank",
      genes = rownames(data$X)
    )
  ))

  # Create desurv_data from newdata
  new_data <- desurv_data(
    X = data$X[, 1:5, drop = FALSE],
    y = data$y[1:5],
    d = data$d[1:5],
    k = 2
  )

  # Should warn about bypassing preprocessing
  expect_warning(
    predict(fit, newdata = new_data),
    "bypasses preprocessing transforms"
  )
})

test_that("predict warns when desurv_data bypasses gene filtering", {
  data <- make_valid_preprocess_data()

  # Create a fit with only gene filtering (no transform)
  fit <- suppressWarnings(desurv_fit(
    X = data$X, y = data$y, d = data$d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 30, ninit = 1,
    verbose = FALSE,
    preprocess_info = list(
      method_trans_train = "none",
      genes = rownames(data$X)[1:5]  # Only a subset of genes
    )
  ))

  # Create desurv_data
  new_data <- desurv_data(
    X = data$X[, 1:5, drop = FALSE],
    y = data$y[1:5],
    d = data$d[1:5],
    k = 2
  )

  # Should warn about bypassing preprocessing (gene filtering)
  expect_warning(
    predict(fit, newdata = new_data),
    "bypasses preprocessing transforms"
  )
})

test_that("predict does not warn for desurv_data when no preprocessing", {
  data <- make_valid_preprocess_data()

  # Create a fit WITHOUT preprocessing metadata
  fit <- suppressWarnings(desurv_fit(
    X = data$X, y = data$y, d = data$d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 30, ninit = 1,
    verbose = FALSE,
    preprocess_info = NULL
  ))

  new_data <- desurv_data(
    X = data$X[, 1:5, drop = FALSE],
    y = data$y[1:5],
    d = data$d[1:5],
    k = 2
  )

  # Should NOT warn when no preprocessing was used
  # Note: there's an existing warning in the code for desurv_data, but only
  # when preprocess_meta has transform or genes. If NULL, no warning.
  expect_no_warning(
    predict(fit, newdata = new_data)
  )
})

# =============================================================================
# Additional validation tests for robustness
# =============================================================================

test_that("preprocess_data() accepts valid data without errors", {
  data <- make_valid_preprocess_data()

  expect_no_error({
    result <- preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    )
  })

  expect_true(!is.null(result$ex))
  expect_true(!is.null(result$sampInfo))
  expect_true(!is.null(result$featInfo))
})

test_that("preprocess_data() works with all-events and all-censored d values", {
  data <- make_valid_preprocess_data()

  # All events
  data$d <- rep(1, length(data$d))
  expect_no_error({
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    )
  })

  # All censored
  data$d <- rep(0, length(data$d))
  expect_no_error({
    preprocess_data(
      X = data$X, y = data$y, d = data$d, dataset = data$dataset,
      method_trans_train = "none", verbose = FALSE
    )
  })
})

test_that("preprocess_data() with rank transform produces expected output", {
  data <- make_valid_preprocess_data()

  result <- preprocess_data(
    X = data$X, y = data$y, d = data$d, dataset = data$dataset,
    method_trans_train = "rank", verbose = FALSE
  )

  # Rank transform should produce values between 1 and nrow(X)
  expect_true(all(result$ex >= 1))
  expect_true(all(result$ex <= nrow(result$ex)))
})

test_that("predict handles non-finite values in newdata for preprocessing path", {
  data <- make_valid_preprocess_data()

  fit <- suppressWarnings(desurv_fit(
    X = data$X, y = data$y, d = data$d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 30, ninit = 1,
    verbose = FALSE,
    preprocess_info = list(
      method_trans_train = "rank",
      genes = rownames(data$X)
    )
  ))

  # Create newdata with Inf - should error at validation
  newX <- data$X[, 1:5, drop = FALSE]
  newX[1, 1] <- Inf

  expect_error(
    predict(fit, newdata = newX),
    "cannot contain NA, NaN, or Inf"
  )
})
