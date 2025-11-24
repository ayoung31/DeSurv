test_that("desurv_bo_refine_bounds respects top-k envelope and importance", {
  bounds <- list(
    k_grid = list(lower = 2, upper = 10, type = "integer"),
    alpha_grid = list(lower = 0, upper = 1)
  )
  history <- data.frame(
    k_grid = c(4, 6, 8, 10),
    alpha_grid = c(0.1, 0.3, 0.6, 0.8),
    mean_cindex = c(0.55, 0.6, 0.58, 0.4),
    status = rep("ok", 4)
  )
  refined <- desurv_bo_refine_bounds(
    bounds,
    history,
    top_k = 2,
    importance = c(k_grid = 1, alpha_grid = 0)
  )
  expect_true(refined$k_grid$upper <= 8)
  expect_true(refined$k_grid$lower >= 4)
  expect_true(refined$alpha_grid$upper <= 0.6)
  expect_true(refined$alpha_grid$lower >= 0.1)
})

test_that("desurv_cv_bayesopt_refine runs with small fixture", {
  skip_on_cran()
  skip_if_not_installed("DiceKriging")
  skip_if_not_installed("lhs")

  fixture <- make_fixture_dataset(n = 16, seed = 77)
  bounds <- list(
    k_grid = list(lower = 2, upper = 3, type = "integer"),
    alpha_grid = list(lower = 0.1, upper = 0.4),
    lambda_grid = list(lower = 0.01, upper = 0.2),
    nu_grid = list(lower = 0.3, upper = 0.7)
  )

  res <- desurv_cv_bayesopt_refine(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    preprocess = FALSE,
    engine = "cold",
    coarse_bounds = bounds,
    bo_fixed = list(lambdaW_grid = 0, lambdaH_grid = 0, nfolds = 2, tol = 1e-4, maxit = 30),
    max_refinements = 1,
    coarse_control = list(n_init = 2, n_iter = 0),
    refine_control = list(n_init = 1, n_iter = 0),
    verbose = FALSE
  )

  expect_s3_class(res, "desurv_cv_bo_refine")
  expect_equal(length(res$runs), 2)
  expect_true(nrow(res$history) >= 2)

  diag <- desurv_bo_refine_diagnostics(res)
  expect_true(nrow(diag$summary) == 2)
})
