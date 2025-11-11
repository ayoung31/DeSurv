test_that("desurv_fit produces coherent factors and predictions", {
  fixture <- make_fixture_dataset()
  fit <- desurv_fit(
    fixture$X,
    fixture$y,
    fixture$d,
    k = fixture$k,
    alpha   = 0.4,
    lambda  = 0.5,
    nu      = 0.5,
    lambdaW = 0.1,
    lambdaH = 0.1,
    ninit   = 2,
    imaxit  = 8,
    maxit   = 60,
    tol     = 1e-4,
    verbose = FALSE
  )

  expect_s3_class(fit, "desurv_fit")
  expect_equal(dim(fit$W), c(fixture$p, fixture$k))
  expect_equal(dim(fit$H), c(fixture$k, fixture$n))
  expect_length(fit$beta, fixture$k)
  expect_true(is.finite(fit$cindex))
  expect_true(fit$cindex >= 0 && fit$cindex <= 1)

  risks <- predict(fit, newdata = fixture$X, type = "risk")
  expect_length(risks, fixture$n)
  expect_true(all(is.finite(risks)))
})

test_that("desurv_alpha_warmstart returns warmstart paths", {
  fixture <- make_fixture_dataset()
  dd <- desurv_data(fixture$X, fixture$y, fixture$d, fixture$k)

  ws <- desurv_alpha_warmstart(
    dd,
    alpha_grid = c(0.2, 0.6),
    lambda     = 0.3,
    nu         = 0.5,
    lambdaW    = 0.05,
    lambdaH    = 0.05,
    n_starts   = 2,
    tol        = 1e-4,
    maxit      = 50,
    verbose    = FALSE
  )

  expect_equal(nrow(ws$results), 4)
  expect_equal(length(ws$fits), 2)
  expect_s3_class(ws$fits[[1]][[1]], "desurv_fit")
})

test_that("desurv_cv selects hyperparameters and refits", {
  fixture <- make_fixture_dataset(n = 36)
  cv_obj <- desurv_cv(
    fixture$X,
    fixture$y,
    fixture$d,
    k_grid       = fixture$k,
    alpha_grid   = c(0.25, 0.5),
    lambda_grid  = 0.4,
    nu_grid      = 0.6,
    lambdaW_grid = 0.05,
    lambdaH_grid = 0.05,
    n_starts     = 1,
    nfolds       = 3,
    seed         = 99,
    tol          = 1e-4,
    maxit        = 80,
    verbose      = FALSE
  )

  expect_s3_class(cv_obj$fit, "desurv_fit")
  expect_true(nrow(cv_obj$summary) > 0)
  expect_setequal(
    names(cv_obj$best),
    c("k", "lambda", "nu", "lambdaW", "lambdaH", "alpha", "mean_cindex", "se_cindex")
  )
  expect_true("cindex" %in% names(cv_obj$fit))
})
