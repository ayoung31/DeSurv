test_that("desurv_bayesopt runs on a small fixture", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 99)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7)
  )

  fit <- desurv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    bo_bounds = bounds,
    bo_fixed = list(
      lambdaW_grid = 0,
      lambdaH_grid = 0,
      nfolds = 2,
      n_starts = 1,
      tol = 1e-4,
      maxit = 30
    ),
    n_init = 2,
    n_iter = 1,
    candidate_pool = 10,
    ntop = 5,
    seed = 123,
    verbose = FALSE
  )

  expect_s3_class(fit, "desurv_bo")
  expect_true(nrow(fit$history) >= 2)
  expect_true(all(names(fit$best$params) %in% names(bounds)))
  expect_false(is.na(fit$best$mean_cindex))
  expect_equal(fit$fixed$ntop, 5)
})

test_that("desurv_bayesopt can warm start from prior history", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 101)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    lambda_grid = list(lower = 0.05, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7)
  )
  bo_fixed <- list(
    lambdaW_grid = 0,
    lambdaH_grid = 0,
    nfolds = 2,
    n_starts = 1,
    tol = 1e-4,
    maxit = 30
  )

  first <- desurv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    bo_bounds = bounds,
    bo_fixed = bo_fixed,
    n_init = 2,
    n_iter = 0,
    candidate_pool = 10,
    ntop = 5,
    seed = 321,
    verbose = FALSE
  )

  second <- desurv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    bo_bounds = bounds,
    bo_fixed = bo_fixed,
    bo_history = first,
    n_init = 1,
    n_iter = 1,
    candidate_pool = 10,
    ntop = 5,
    seed = 654,
    verbose = FALSE
  )

  expect_true(nrow(second$history) > nrow(first$history))
  expect_equal(
    second$history[seq_len(nrow(first$history)), names(bounds)],
    first$history[, names(bounds)]
  )
  expect_true(any(second$history$stage == "bo"))
})

test_that("desurv_bayesopt warm start handles factor histories", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 55)

  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.1, upper = 0.3),
    lambda_grid = list(lower = 0.01, upper = 0.2),
    nu_grid = list(lower = 0.4, upper = 0.6)
  )
  bo_fixed <- list(
    lambdaW_grid = 0,
    lambdaH_grid = 0,
    nfolds = 2,
    n_starts = 1,
    tol = 1e-4,
    maxit = 30
  )

  first <- desurv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    bo_bounds = bounds,
    bo_fixed = bo_fixed,
    n_init = 2,
    n_iter = 0,
    candidate_pool = 10,
    ntop = 5,
    seed = 999,
    verbose = FALSE
  )

  hist_df <- first$history
  num_cols <- names(bounds)
  for (nm in num_cols) {
    hist_df[[nm]] <- factor(hist_df[[nm]])
  }
  hist_df$stage <- factor(hist_df$stage)
  hist_df$iteration <- factor(hist_df$iteration)
  hist_df$mean_cindex <- factor(sprintf("%.6f", hist_df$mean_cindex))
  hist_df$status <- factor(hist_df$status)
  hist_df$message <- factor(hist_df$message)
  hist_df$elapsed <- factor(hist_df$elapsed)

  warm <- desurv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    bo_bounds = bounds,
    bo_fixed = bo_fixed,
    bo_history = hist_df,
    n_init = 1,
    n_iter = 1,
    candidate_pool = 10,
    ntop = 5,
    seed = 998,
    verbose = FALSE
  )

  expect_gt(nrow(warm$history), nrow(first$history))
  expect_true(any(warm$history$stage == "bo"))
})

test_that("desurv_bo_default_bounds exposes expected keys", {
  bounds <- desurv_bo_default_bounds(include_factor_penalties = FALSE)
  expect_type(bounds, "list")
  expect_true(all(c("k_grid", "alpha_grid", "lambda_grid") %in% names(bounds)))
  expect_false(any(grepl("lambdaW_grid", names(bounds))))
  expect_true("ntop" %in% names(bounds))
})

# Tests for Finding 5: desurv_bo_default_bounds() include_ngene parameter
test_that("desurv_bo_default_bounds include_ngene parameter controls ngene inclusion", {
  # With include_ngene = TRUE (default)
  bounds_with_ngene <- desurv_bo_default_bounds(include_ngene = TRUE)
  expect_true("ngene" %in% names(bounds_with_ngene))

  # With include_ngene = FALSE
  bounds_without_ngene <- desurv_bo_default_bounds(include_ngene = FALSE)
  expect_false("ngene" %in% names(bounds_without_ngene))

  # Should still have other core parameters
  expect_true("k_grid" %in% names(bounds_without_ngene))
  expect_true("alpha_grid" %in% names(bounds_without_ngene))
})

# Tests for Finding 2: Fixed parameters exclude tuned parameters
test_that("desurv_bayesopt $fixed excludes parameters in bo_bounds", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 42)

  # Include nfolds in bo_bounds - it should NOT appear in $fixed
  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.4),
    nfolds = list(lower = 2, upper = 3, type = "integer")  # Tuned, not fixed
  )

  fit <- suppressWarnings(desurv_bayesopt(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    bo_bounds = bounds,
    bo_fixed = list(
      lambdaW_grid = 0,
      lambdaH_grid = 0,
      lambda_grid = 0.1,
      nu_grid = 0.5,
      n_starts = 1,
      tol = 1e-4,
      maxit = 30
    ),
    n_init = 2,
    n_iter = 0,
    candidate_pool = 5,
    seed = 123,
    verbose = FALSE
  ))

  # nfolds is in bo_bounds, so should NOT be in $fixed
  expect_false("nfolds" %in% names(fit$fixed))

  # lambdaW_grid is in bo_fixed but not bo_bounds, so should be in $fixed
  expect_true("lambdaW_grid" %in% names(fit$fixed))
})

# Tests for Finding 4: desurv_bo_refine_bounds validates top_k
test_that("desurv_bo_refine_bounds validates top_k as positive integer", {
  bounds <- list(
    k_grid = list(lower = 2, upper = 5, type = "integer"),
    alpha_grid = list(lower = 0.1, upper = 0.5)
  )
  history <- data.frame(
    k_grid = c(2, 3, 4),
    alpha_grid = c(0.2, 0.3, 0.4),
    mean_cindex = c(0.6, 0.7, 0.65),
    status = c("ok", "ok", "ok"),
    stringsAsFactors = FALSE
  )

  # top_k = 0 should error
  expect_error(
    desurv_bo_refine_bounds(bounds, history, top_k = 0),
    "positive integer"
  )

  # top_k = -1 should error
  expect_error(
    desurv_bo_refine_bounds(bounds, history, top_k = -1),
    "positive integer"
  )

  # top_k = NA should error
  expect_error(
    desurv_bo_refine_bounds(bounds, history, top_k = NA),
    "positive integer"
  )

  # top_k = 1 should work
  result <- desurv_bo_refine_bounds(bounds, history, top_k = 1)
  expect_type(result, "list")
})

# Tests for Finding 3: Failed evaluations don't block regions in warm-start
test_that("warm-start history only blocks successful evaluation keys", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  # Create a mock history with one successful and one failed evaluation
  # The failed one should NOT block that region in warm-start
  history_df <- data.frame(
    eval_id = c(1, 2),
    stage = c("initial", "initial"),
    iteration = c(0, 0),
    k_grid = c(2, 3),
    alpha_grid = c(0.25, 0.35),
    mean_cindex = c(0.65, NA),  # Second one failed
    status = c("ok", "error"),
    message = c(NA, "test error"),
    elapsed = c(1.0, 0.5),
    stringsAsFactors = FALSE
  )

  bounds <- list(
    k_grid = list(lower = 2, upper = 4, type = "integer"),
    alpha_grid = list(lower = 0.2, upper = 0.5)
  )

  bound_info <- DeSurv:::.desurv_bo_parse_bounds(bounds)
  warmstart <- DeSurv:::.desurv_bo_prepare_history(history_df, bound_info)

  # Should have imported both rows (for history display)
  expect_equal(warmstart$n_loaded, 2L)

  # But existing_keys should only block the successful one
  # The failed one (k_grid=3, alpha_grid=0.35) should be retryable
  expect_equal(length(warmstart$existing_keys), 1L)
})

# Tests for desurv_cv_bo_default_bounds parameters
test_that("desurv_cv_bo_default_bounds respects include_ngene parameter", {
  # With include_ngene = TRUE (default)
  bounds_with <- desurv_cv_bo_default_bounds(include_ngene = TRUE)
  expect_true("ngene" %in% names(bounds_with))

  # With include_ngene = FALSE
  bounds_without <- desurv_cv_bo_default_bounds(include_ngene = FALSE)
  expect_false("ngene" %in% names(bounds_without))
})
