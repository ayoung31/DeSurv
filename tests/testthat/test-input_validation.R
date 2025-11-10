test_that("input_validation passes on valid input", {
  set.seed(1)
  p <- 10; n <- 20; k <- 3
  X <- matrix(rexp(p * n), nrow = p, ncol = n)
  y <- rexp(n)
  d <- rbinom(n, 1, 0.6)
  W0 <- matrix(abs(rnorm(p * k)), nrow = p)
  H0 <- matrix(abs(rnorm(k * n)), nrow = k)
  beta0 <- rnorm(k)

  res <- input_validation(
    X, y, d, k, W0, H0, beta0,
    alpha = 0.5, lambda = 1, nu = 0.3,
    lambdaW = 0.1, lambdaH = 0,
    tol = 1e-5, maxit = 200, verbose = TRUE
  )

  expect_type(res, "list")
  expect_equal(res$k, k)
  expect_true(all(sapply(res[c("alpha","nu","lambda","lambdaW","lambdaH","tol","maxit")], is.numeric)))
})

test_that("X must be numeric, finite, and non-negative", {
  p <- 5; n <- 6; k <- 2
  X <- matrix(abs(rnorm(p * n)), nrow = p)
  y <- rexp(n); d <- rbinom(n, 1, 0.5)
  W0 <- matrix(abs(rnorm(p * k)), p)
  H0 <- matrix(abs(rnorm(k * n)), k)
  beta0 <- rnorm(k)

  expect_error(input_validation(as.list(X), y, d, k, W0, H0, beta0,
                                0.5, 1, 0.5, 0.1, 0.1, 1e-4, 100, FALSE),
               "numeric matrix")
  expect_error(input_validation(-X, y, d, k, W0, H0, beta0,
                                0.5, 1, 0.5, 0.1, 0.1, 1e-4, 100, FALSE),
               "non-negative")
  X_bad <- X; X_bad[1,1] <- NA
  expect_error(input_validation(X_bad, y, d, k, W0, H0, beta0,
                                0.5, 1, 0.5, 0.1, 0.1, 1e-4, 100, FALSE),
               "NA, NaN, or Inf")
})

test_that("dimension and survival checks hold", {
  p <- 5; n <- 6; k <- 2
  X <- matrix(abs(rnorm(p * n)), p)
  y <- rexp(n); d <- rbinom(n, 1, 0.5)
  W0 <- matrix(abs(rnorm(p * k)), p)
  H0 <- matrix(abs(rnorm(k * n)), k)
  beta0 <- rnorm(k)

  # mismatched dimensions
  expect_error(input_validation(X, y, d, k, W0[-1,,drop=FALSE], H0, beta0,
                                0.5,1,0.5,0.1,0.1,1e-4,100,FALSE),
               "nrow\\(W0\\)")
  expect_error(input_validation(X, y, d, k, W0, H0[, -1, drop=FALSE], beta0,
                                0.5,1,0.5,0.1,0.1,1e-4,100,FALSE),
               "ncol\\(H0\\)")

  # survival vector issues
  expect_error(input_validation(X, y[-1], d, k, W0, H0, beta0,
                                0.5,1,0.5,0.1,0.1,1e-4,100,FALSE),
               "match number of columns")
  expect_error(input_validation(X, y, rep(0,n), k, W0, H0, beta0,
                                0.5,1,0.5,0.1,0.1,1e-4,100,FALSE),
               "must be > 0")
})

test_that("hyperparameters clip and coerce correctly", {
  p <- 5; n <- 6; k <- 2
  X <- matrix(abs(rnorm(p * n)), p)
  y <- rexp(n); d <- rbinom(n, 1, 0.5)
  W0 <- matrix(abs(rnorm(p * k)), p)
  H0 <- matrix(abs(rnorm(k * n)), k)
  beta0 <- rnorm(k)

  expect_warning(res <- input_validation(X, y, d, k, W0, H0, beta0,
                                         1.5, 1, 0.5, 0.1, 0.1,
                                         1e-4, 100, FALSE),
                 "alpha outside")
  expect_lt(res$alpha, 1)

  expect_warning(res <- input_validation(X, y, d, k, W0, H0, beta0,
                                         0.5, 1, -0.5, 0.1, 0.1,
                                         1e-4, 100, FALSE),
                 "nu outside")
  expect_gte(res$nu, 0)
  expect_lte(res$nu, 1)

  expect_warning(res <- input_validation(X, y, d, k, W0, H0, beta0,
                                         0.5, -1, 0.5, -0.1, -0.2,
                                         1e-4, 100, FALSE),
                 "negative")
  expect_equal(res$lambda, 0)
  expect_equal(res$lambdaW, 0)
  expect_equal(res$lambdaH, 0)
})
