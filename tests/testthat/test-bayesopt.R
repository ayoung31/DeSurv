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

test_that("desurv_bo_default_bounds exposes expected keys", {
  bounds <- desurv_bo_default_bounds(include_factor_penalties = FALSE)
  expect_type(bounds, "list")
  expect_true(all(c("k_grid", "alpha_grid", "lambda_grid") %in% names(bounds)))
  expect_false(any(grepl("lambdaW_grid", names(bounds))))
  expect_true("ntop" %in% names(bounds))
})
