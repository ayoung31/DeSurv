test_that("desurv_cv_bayesopt runs on a small fixture", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 99)

  # Define a small search space for testing
  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7),
    lambdaW_grid = list(lower = 0, upper = 0.01),
    lambdaH_grid = list(lower = 0, upper = 0.01),
    n_starts = list(lower = 1, upper = 2, type = "integer"),
    ntop = list(lower = 10, upper = 20, type = "integer")
  )

  fit <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 2,
    tol = 1e-4,
    maxit = 30,
    bo_bounds = bounds,
    n_init = 2,
    n_iter = 1,
    candidate_pool = 10,
    seed = 123,
    verbose = FALSE
  )

  expect_s3_class(fit, "desurv_cv_bo")
  expect_true(nrow(fit$history) >= 2)
  expect_true(all(names(fit$best$params) %in% names(bounds)))
  expect_false(is.na(fit$best$mean_cindex))
  expect_true(fit$best$mean_cindex >= 0 && fit$best$mean_cindex <= 1)
  expect_true("ntop" %in% names(fit$history))
  expect_null(fit$fixed$ntop)
})

test_that("desurv_cv_bayesopt can warm start from a data frame", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 109)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7),
    lambdaW_grid = list(lower = 0, upper = 0.01),
    lambdaH_grid = list(lower = 0, upper = 0.01),
    n_starts = list(lower = 1, upper = 2, type = "integer"),
    ntop = list(lower = 10, upper = 20, type = "integer")
  )

  first <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 2,
    tol = 1e-4,
    maxit = 30,
    bo_bounds = bounds,
    n_init = 2,
    n_iter = 0,
    candidate_pool = 10,
    seed = 111,
    verbose = FALSE
  )

  second <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 2,
    tol = 1e-4,
    maxit = 30,
    bo_bounds = bounds,
    bo_history = first$history,
    n_init = 1,
    n_iter = 1,
    candidate_pool = 10,
    seed = 222,
    verbose = FALSE
  )

  expect_true(nrow(second$history) > nrow(first$history))
  expect_equal(
    second$history[seq_len(nrow(first$history)), names(bounds)],
    first$history[, names(bounds)]
  )
  expect_true(any(second$history$stage == "bo"))
})

test_that("desurv_cv_bayesopt warm start handles factor histories", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 404)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.1, upper = 0.3),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.4, upper = 0.6),
    lambdaW_grid = list(lower = 0, upper = 0.01),
    lambdaH_grid = list(lower = 0, upper = 0.01)
  )

  first <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 2,
    tol = 1e-4,
    maxit = 30,
    bo_bounds = bounds,
    n_init = 2,
    n_iter = 0,
    candidate_pool = 10,
    seed = 1111,
    verbose = FALSE
  )

  hist_df <- first$history
  for (nm in names(bounds)) {
    hist_df[[nm]] <- factor(hist_df[[nm]])
  }
  hist_df$stage <- factor(hist_df$stage)
  hist_df$iteration <- factor(hist_df$iteration)
  hist_df$mean_cindex <- factor(sprintf("%.6f", hist_df$mean_cindex))
  hist_df$status <- factor(hist_df$status)
  hist_df$message <- factor(hist_df$message)
  hist_df$elapsed <- factor(hist_df$elapsed)

  warm <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 2,
    tol = 1e-4,
    maxit = 30,
    bo_bounds = bounds,
    bo_history = hist_df,
    n_init = 1,
    n_iter = 1,
    candidate_pool = 10,
    seed = 2222,
    verbose = FALSE
  )

  expect_gt(nrow(warm$history), nrow(first$history))
  expect_true(any(warm$history$stage == "bo"))
})

test_that("desurv_cv_bo_default_bounds exposes expected keys", {
  bounds <- desurv_cv_bo_default_bounds(
    include_factor_penalties = FALSE,
    include_n_starts = TRUE,
    include_ngene = TRUE
  )
  expect_type(bounds, "list")
  expect_true(all(c("k_grid", "alpha_grid", "lambda_grid", "nu_grid") %in% names(bounds)))
  expect_false(any(grepl("lambdaW_grid", names(bounds))))
  expect_true("n_starts" %in% names(bounds))
  expect_true("ngene" %in% names(bounds))
  expect_true("ntop" %in% names(bounds))
})

test_that("desurv_cv_bo_default_bounds respects include flags", {
  # Test with all includes FALSE except core params
  bounds_minimal <- desurv_cv_bo_default_bounds(
    include_factor_penalties = FALSE,
    include_n_starts = FALSE,
    include_ngene = FALSE
  )

  expect_true(all(c("k_grid", "alpha_grid", "lambda_grid", "nu_grid") %in% names(bounds_minimal)))
  expect_false("lambdaW_grid" %in% names(bounds_minimal))
  expect_false("lambdaH_grid" %in% names(bounds_minimal))
  expect_false("n_starts" %in% names(bounds_minimal))
  expect_false("ngene" %in% names(bounds_minimal))
  expect_true("ntop" %in% names(bounds_minimal))

  # Test with all includes TRUE
  bounds_full <- desurv_cv_bo_default_bounds(
    include_factor_penalties = TRUE,
    include_n_starts = TRUE,
    include_ngene = TRUE
  )

  expect_true("lambdaW_grid" %in% names(bounds_full))
  expect_true("lambdaH_grid" %in% names(bounds_full))
  expect_true("n_starts" %in% names(bounds_full))
  expect_true("ngene" %in% names(bounds_full))
  expect_true("ntop" %in% names(bounds_full))
})

test_that("desurv_cv_bayesopt handles errors gracefully", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 99)

  # Empty bounds should error
  expect_error(
    desurv_cv_bayesopt(
      X = fixture$X,
      y = fixture$y,
      d = fixture$d,
      bo_bounds = list(),
      verbose = FALSE
    ),
    "must contain at least one tunable parameter"
  )
})

test_that("desurv_cv_bayesopt fixed params are stored correctly", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 99)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7),
    lambdaW_grid = list(lower = 0, upper = 0.01),
    lambdaH_grid = list(lower = 0, upper = 0.01),
    n_starts = list(lower = 1, upper = 2, type = "integer"),
    ntop = list(lower = 10, upper = 20, type = "integer")
  )

  fit <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 3,
    tol = 1e-4,
    maxit = 100,
    bo_bounds = bounds,
    n_init = 2,
    n_iter = 1,
    seed = 123,
    verbose = FALSE
  )

  expect_equal(fit$fixed$nfolds, 3)
  expect_equal(fit$fixed$tol, 1e-4)
  expect_equal(fit$fixed$maxit, 100)
  expect_equal(fit$fixed$engine, "cold")
  expect_null(fit$fixed$ntop)
})

test_that("print.desurv_cv_bo works", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 99)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7),
    lambdaW_grid = list(lower = 0, upper = 0.01),
    lambdaH_grid = list(lower = 0, upper = 0.01),
    n_starts = list(lower = 1, upper = 2, type = "integer"),
    ntop = list(lower = 10, upper = 20, type = "integer")
  )

  fit <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 2,
    tol = 1e-4,
    maxit = 100,
    bo_bounds = bounds,
    n_init = 2,
    n_iter = 1,
    seed = 123,
    verbose = FALSE
  )

  # Should not error
  expect_output(print(fit), "Bayesian optimisation for desurv_cv")
  expect_output(print(fit), "Evaluations")
  expect_output(print(fit), "Best C-index")
})

test_that("desurv_cv_bayesopt stores ntop when provided", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 99)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7),
    lambdaW_grid = list(lower = 0, upper = 0.01),
    lambdaH_grid = list(lower = 0, upper = 0.01),
    n_starts = list(lower = 1, upper = 2, type = "integer")
  )

  fit <- desurv_cv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    nfolds = 2,
    tol = 1e-4,
    maxit = 30,
    ntop = 40,
    bo_bounds = bounds,
    n_init = 2,
    n_iter = 1,
    seed = 456,
    verbose = FALSE
  )

  expect_equal(fit$fixed$ntop, 40)
})
