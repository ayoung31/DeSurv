test_that("desurv_fit produces coherent factors and predictions", {
  fixture <- make_fixture_dataset()
  fit <- suppressWarnings(
    desurv_fit(
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

  ws <- suppressWarnings(
    desurv_alpha_warmstart(
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
  )

  expect_equal(nrow(ws$results), 4)
  expect_equal(length(ws$fits), 2)
  expect_s3_class(ws$fits[[1]][[1]], "desurv_fit")
})

test_that("desurv_cv selects hyperparameters and refits", {
  fixture <- make_fixture_dataset(n = 36)
  cv_obj <- suppressWarnings(
    desurv_cv(
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
  )

  expect_s3_class(cv_obj$fit, "desurv_fit")
  expect_true(nrow(cv_obj$summary) > 0)
  expect_setequal(
    names(cv_obj$best),
    c("k", "lambda", "nu", "lambdaW", "lambdaH", "alpha", "mean_cindex", "se_cindex")
  )
  expect_true("cindex" %in% names(cv_obj$fit))
})

test_that("desurv_cv cold engine returns stable CV summaries", {
  fixture <- make_fixture_dataset(n = 32, seed = 88)
  alpha_grid <- c(0.2, 0.5)

  base_args <- list(
    X            = fixture$X,
    y            = fixture$y,
    d            = fixture$d,
    k_grid       = fixture$k,
    alpha_grid   = alpha_grid,
    lambda_grid  = 0.5,
    nu_grid      = 0.4,
    lambdaW_grid = 0.05,
    lambdaH_grid = 0.05,
    n_starts     = 2,
    nfolds       = 3,
    tol          = 1e-4,
    maxit        = 150,
    verbose      = FALSE
  )

  cold_1 <- suppressWarnings(do.call(
    desurv_cv,
    c(base_args, list(seed = 2024, engine = "cold"))
  ))
  cold_2 <- suppressWarnings(do.call(
    desurv_cv,
    c(base_args, list(seed = 2024, engine = "cold"))
  ))

  expect_s3_class(cold_1$fit, "desurv_fit")
  expect_equal(cold_1$best, cold_2$best)
  expect_equal(cold_1$summary, cold_2$summary)

  expected_rows <- length(alpha_grid) *
    length(base_args$k_grid) *
    length(base_args$lambda_grid) *
    length(base_args$nu_grid) *
    length(base_args$lambdaW_grid) *
    length(base_args$lambdaH_grid) *
    base_args$nfolds *
    base_args$n_starts
  expect_equal(nrow(cold_1$cv_results), expected_rows)
  expect_true(all(c("alpha", "cindex") %in% names(cold_1$cv_results)))
  expect_equal(length(unique(cold_1$cv_results$init_id)), base_args$n_starts)
})

test_that("desurv_cv cv_only returns raw CV object without refitting", {
  fixture <- make_fixture_dataset(n = 24, seed = 22)

  cv_obj <- suppressWarnings(
    desurv_cv(
      fixture$X,
      fixture$y,
      fixture$d,
      k_grid       = fixture$k,
      alpha_grid   = c(0.2, 0.4),
      lambda_grid  = 0.3,
      nu_grid      = 0.5,
      lambdaW_grid = 0.05,
      lambdaH_grid = 0.05,
      n_starts     = 1,
      nfolds       = 2,
      tol          = 1e-4,
      maxit        = 60,
      verbose      = FALSE,
      cv_only      = TRUE
    )
  )

  expect_s3_class(cv_obj, "desurv_cv")
  expect_true(is.null(cv_obj$fit))
  expect_true(is.data.frame(cv_obj$summary))
  expect_equal(cv_obj$rule, "max")
})

test_that("desurv_cv warns and reduces folds when nfolds exceeds data capacity", {
  fixture <- make_fixture_dataset(n = 12, seed = 10)

  expect_warning(
    fit <- suppressWarnings(
      desurv_cv(
        fixture$X,
        fixture$y,
        fixture$d,
        k_grid       = fixture$k,
        alpha_grid   = 0.4,
        lambda_grid  = 0.3,
        nu_grid      = 0.5,
        lambdaW_grid = 0.05,
        lambdaH_grid = 0.05,
        n_starts     = 1,
        nfolds       = 20,
        tol          = 1e-4,
        maxit        = 50,
        verbose      = FALSE
      )
    ),
    "reducing nfolds"
  )

  folds_used <- sort(unique(fit$folds))
  expect_equal(folds_used, seq_along(folds_used))
  expect_equal(
    sort(unique(fit$cv_results$fold)),
    folds_used
  )
})

test_that("summary retains folds whose C-index is undefined", {
  X <- matrix(rexp(4 * 4), nrow = 4)
  y <- c(5, 4, 3, 2)
  d <- c(1, 1, 0, 0)
  folds <- c(1, 1, 2, 2)

  fit <- suppressWarnings(
    desurv_cv(
      X,
      y,
      d,
      k_grid       = 2,
      alpha_grid   = 0.3,
      lambda_grid  = 0.4,
      nu_grid      = 0.5,
      lambdaW_grid = 0.05,
      lambdaH_grid = 0.05,
      n_starts     = 1,
      nfolds       = 2,
      folds        = folds,
      tol          = 1e-4,
      maxit        = 40,
      verbose      = FALSE
    )
  )

  fold2 <- subset(fit$summary_fold, fold == 2)
  expect_equal(nrow(fold2), length(unique(fit$summary_fold$alpha)))
  expect_true(all(is.na(fold2$mean_cindex_fold)))
  expect_true(any(is.na(fit$summary$mean_cindex)))
})

test_that("predict.desurv_fit enforces gene alignment and rowname requirements", {
  fixture <- make_fixture_dataset()
  fit <- suppressWarnings(
    desurv_fit(
      fixture$X,
      fixture$y,
      fixture$d,
      k = fixture$k,
      alpha   = 0.3,
      lambda  = 0.4,
      nu      = 0.5,
      lambdaW = 0.1,
      lambdaH = 0.1,
      ninit   = 2,
      imaxit  = 5,
      maxit   = 40,
      tol     = 1e-4,
      verbose = FALSE
    )
  )

  shuffled <- fixture$X[rev(rownames(fixture$X)), , drop = FALSE]
  preds_ref <- predict(fit, type = "lp")
  preds_perm <- predict(fit, newdata = shuffled, type = "lp")
  expect_equal(preds_perm, preds_ref, tolerance = 1e-12)

  missing_gene <- shuffled[-1, , drop = FALSE]
  expect_error(
    predict(fit, newdata = missing_gene),
    "missing",
    fixed = FALSE
  )

  unlabelled <- shuffled
  rownames(unlabelled) <- NULL
  expect_error(
    predict(fit, newdata = unlabelled),
    "rownames",
    fixed = FALSE
  )
})
