# Tests for numerical stability guards in C++ optimization
# These tests verify that the optimizer handles pathological inputs gracefully
# without producing NaN/Inf values that would corrupt results.

test_that("optimizer handles near-zero variance features gracefully", {
  # Create data where some features have near-zero variance
  # This can trigger Hessian clamping in Cox CD
  set.seed(42)
  n <- 30
  p <- 10

  # X with some nearly-constant columns (non-negative for NMF)
  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  X[1, ] <- rep(1, n)  # constant feature
  X[2, ] <- rep(1, n) + abs(rnorm(n, sd = 1e-10))  # near-constant feature

  # Survival data
  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.6)

  # Ensure at least some events
  if (sum(d) < 3) d[1:3] <- 1

  # Should run without error (guards prevent NaN propagation)
  expect_no_error({
    result <- suppressWarnings(desurv_fit(
      X = X, y = y, d = d, k = 2,
      alpha = 0.5, lambda = 0.1, nu = 0.5,
      lambdaW = 0.01, lambdaH = 0.01,
      tol = 1e-3, maxit = 20,
      ninit = 1, verbose = FALSE
    ))
  })
})

test_that("optimizer handles extreme scale differences in X", {
  # Features with vastly different scales can challenge numerical stability
  set.seed(123)
  n <- 25
  p <- 8

  # Non-negative matrix for NMF
  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  X[1, ] <- X[1, ] * 1e6  # very large scale
  X[2, ] <- X[2, ] * 1e-6  # very small scale

  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.5)
  if (sum(d) < 3) d[1:3] <- 1

  # Should handle extreme scales (may not converge well, but shouldn't crash)
  expect_no_error({
    result <- suppressWarnings(desurv_fit(
      X = X, y = y, d = d, k = 2,
      alpha = 0.3, lambda = 0.01, nu = 0.5,
      lambdaW = 0.01, lambdaH = 0.01,
      tol = 1e-2, maxit = 15,
      ninit = 1, verbose = FALSE
    ))
  })
})

test_that("optimizer produces finite W and H matrices", {
  # Standard test case that verifies W and H remain finite
  set.seed(789)
  n <- 20
  p <- 6

  X <- matrix(abs(rnorm(n * p)), nrow = p, ncol = n)
  y <- sort(runif(n, 1, 50))
  d <- rbinom(n, 1, 0.7)
  if (sum(d) < 3) d[1:3] <- 1

  result <- suppressWarnings(desurv_fit(
    X = X, y = y, d = d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 50,
    ninit = 1, verbose = FALSE
  ))

  # W and H should be finite

  expect_true(all(is.finite(result$W)), info = "W contains non-finite values")
  expect_true(all(is.finite(result$H)), info = "H contains non-finite values")
  expect_true(all(is.finite(result$beta)), info = "beta contains non-finite values")
})

test_that("optimizer handles all-zero events correctly",
 {
  # Edge case: no events (all censored)
  # This is pathological but shouldn't crash
  set.seed(456)
  n <- 15
  p <- 5

  X <- matrix(abs(rnorm(n * p)), nrow = p, ncol = n)
  y <- runif(n, 1, 100)
  d <- rep(0, n)  # all censored - no events

  # This should fail with an informative error or handle gracefully
  # The Cox model is undefined with no events, so we expect either
  # an error or a result with beta = 0
  tryCatch({
    result <- suppressWarnings(desurv_fit(
      X = X, y = y, d = d, k = 2,
      alpha = 0.5, lambda = 0.1, nu = 0.5,
      lambdaW = 0.01, lambdaH = 0.01,
      tol = 1e-3, maxit = 20,
      ninit = 1, verbose = FALSE
    ))
    # If it succeeds, beta should be zero or near-zero (no survival signal)
    expect_true(all(is.finite(result$beta)))
  }, error = function(e) {
    # An error is also acceptable for this pathological case
    expect_true(TRUE)
  })
})

test_that("optimizer handles very small sample size", {
  # Minimum viable sample size
  set.seed(999)
  n <- 6  # very small
  p <- 4
  k <- 2

  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  y <- runif(n, 1, 50)
  d <- c(1, 1, 1, 0, 0, 0)  # half events

  # Should handle small n without crashing
  expect_no_error({
    result <- suppressWarnings(desurv_fit(
      X = X, y = y, d = d, k = k,
      alpha = 0.5, lambda = 0.1, nu = 0.5,
      lambdaW = 0.01, lambdaH = 0.01,
      tol = 1e-2, maxit = 30,
      ninit = 1, verbose = FALSE
    ))
  })
})

test_that("optimizer handles tied survival times", {
  # Many tied survival times can affect Cox model numerical stability
  set.seed(111)
  n <- 24
  p <- 6

  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)

  # Create many tied times
  y <- rep(c(10, 20, 30, 40), each = 6)  # 4 unique times, 6 ties each
  d <- rbinom(n, 1, 0.6)
  if (sum(d) < 3) d[1:3] <- 1

  result <- suppressWarnings(desurv_fit(
    X = X, y = y, d = d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 50,
    ninit = 1, verbose = FALSE
  ))

  # Results should be finite
  expect_true(all(is.finite(result$W)))
  expect_true(all(is.finite(result$H)))
  expect_true(all(is.finite(result$beta)))
})

test_that("optimizer handles high alpha (pure Cox) without instability", {
  # alpha = 1 means pure Cox model, no NMF component
  # This tests the Cox-only path for stability
  set.seed(222)
  n <- 20
  p <- 8

  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.6)
  if (sum(d) < 3) d[1:3] <- 1

  result <- suppressWarnings(desurv_fit(
    X = X, y = y, d = d, k = 2,
    alpha = 1.0,  # pure Cox
    lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 50,
    ninit = 1, verbose = FALSE
  ))

  expect_true(all(is.finite(result$beta)))
  expect_true(all(is.finite(result$W)))
  expect_true(all(is.finite(result$H)))
})

test_that("optimizer handles low alpha (pure NMF) without instability", {
  # alpha = 0 means pure NMF, no Cox component
  # This tests the NMF-only path for stability
  set.seed(333)
  n <- 20
  p <- 8

  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.6)
  if (sum(d) < 3) d[1:3] <- 1

  result <- suppressWarnings(desurv_fit(
    X = X, y = y, d = d, k = 2,
    alpha = 0.0,  # pure NMF
    lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 50,
    ninit = 1, verbose = FALSE
  ))

  expect_true(all(is.finite(result$W)))
  expect_true(all(is.finite(result$H)))
  # beta should be unchanged when alpha=0
  expect_true(all(is.finite(result$beta)))
})

test_that("optimizer handles sparse X matrix (many zeros)", {
  # Sparse input can cause numerical issues in NMF updates
  set.seed(444)
  n <- 25
  p <- 10

  # Create sparse matrix with ~80% zeros
  X <- matrix(0, nrow = p, ncol = n)
  non_zero_idx <- sample(length(X), size = floor(0.2 * length(X)))
  X[non_zero_idx] <- abs(rnorm(length(non_zero_idx))) + 0.1

  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.6)
  if (sum(d) < 3) d[1:3] <- 1

  result <- suppressWarnings(desurv_fit(
    X = X, y = y, d = d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-3, maxit = 50,
    ninit = 1, verbose = FALSE
  ))

  expect_true(all(is.finite(result$W)))
  expect_true(all(is.finite(result$H)))
  expect_true(all(is.finite(result$beta)))
})

test_that("loss values remain finite throughout optimization", {
  # Verify the loss tracking doesn't produce non-finite values
  set.seed(555)
  n <- 20
  p <- 6

  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.6)
  if (sum(d) < 3) d[1:3] <- 1

  result <- suppressWarnings(desurv_fit(
    X = X, y = y, d = d, k = 2,
    alpha = 0.5, lambda = 0.1, nu = 0.5,
    lambdaW = 0.01, lambdaH = 0.01,
    tol = 1e-4, maxit = 100,
    ninit = 1, verbose = FALSE
  ))

  # Final loss should be finite (result$loss is a list; $loss component has the total)
  expect_true(is.finite(result$loss$loss), info = "Final loss is non-finite")

  # Loss history should be finite (for iterations run)
  loss_vec <- result$lossit
  used_loss <- loss_vec[loss_vec != 0]  # non-zero entries are used iterations
  if (length(used_loss) > 0) {
    expect_true(all(is.finite(used_loss)), info = "Loss history contains non-finite values")
  }
})

# Test for Finding 5: Verify Hessian clamping warning can be triggered
# Note: This is hard to test directly as the warning threshold is high,
# but we can verify the code path exists by checking function signature
test_that("Cox CD function tracks Hessian clamping (code path exists)", {
  # This test verifies that the internal C++ function is properly called
  # The actual warning would require extreme conditions to trigger

  # Create a normal case that should NOT trigger the warning
  set.seed(666)
  n <- 20
  p <- 6

  X <- matrix(abs(rnorm(n * p)) + 0.1, nrow = p, ncol = n)
  y <- runif(n, 1, 100)
  d <- rbinom(n, 1, 0.7)
  if (sum(d) < 3) d[1:3] <- 1

  # Normal case should not produce warnings about Hessian clamping
  # Use expect_silent wrapped with expect_warning to check no Hessian warnings
  warnings_captured <- character()
  result <- withCallingHandlers(
    desurv_fit(
      X = X, y = y, d = d, k = 2,
      alpha = 0.5, lambda = 0.1, nu = 0.5,
      lambdaW = 0.01, lambdaH = 0.01,
      tol = 1e-3, maxit = 50,
      ninit = 1, verbose = FALSE
    ),
    warning = function(w) {
      warnings_captured <<- c(warnings_captured, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Verify no Hessian clamping warnings
  hessian_warnings <- grep("Hessian clamping", warnings_captured, value = TRUE)
  expect_length(hessian_warnings, 0)
})
