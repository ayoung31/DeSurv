test_that("predict applies rank preprocessing metadata to raw matrices", {
  fixture <- make_fixture_dataset(n = 12, seed = 101)
  dataset <- rep("A", fixture$n)

  prep <- preprocess_data(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    dataset = dataset,
    ngene = fixture$p,
    method_trans_train = "rank",
    verbose = FALSE
  )

  preprocess_meta <- list(
    genes = prep$featInfo,
    method_trans_train = "rank"
  )

  fit <- suppressWarnings(desurv_fit(
    X       = prep$ex,
    y       = prep$sampInfo$time,
    d       = prep$sampInfo$event,
    k       = fixture$k,
    alpha   = 0.3,
    lambda  = 0.2,
    nu      = 0.5,
    lambdaW = 0.05,
    lambdaH = 0.05,
    verbose = FALSE,
    preprocess_info = preprocess_meta
  ))

  baseline <- predict(
    fit,
    newdata = desurv_data(prep$ex, prep$sampInfo$time, prep$sampInfo$event, fixture$k)
  )

  raw_subset <- fixture$X[prep$featInfo, , drop = FALSE]
  pred_raw <- predict(fit, newdata = raw_subset)
  expect_equal(pred_raw, baseline)
})

test_that("predict applies quantile preprocessing metadata to raw matrices", {
  skip_if_not_installed("preprocessCore")
  fixture <- make_fixture_dataset(n = 12, seed = 77)
  dataset <- rep("A", fixture$n)

  prep <- preprocess_data(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    dataset = dataset,
    ngene = fixture$p,
    method_trans_train = "quant",
    verbose = FALSE
  )

  preprocess_meta <- list(
    genes = prep$featInfo,
    method_trans_train = "quant",
    transform_target = prep$transform_target
  )

  fit <- suppressWarnings(desurv_fit(
    X       = prep$ex,
    y       = prep$sampInfo$time,
    d       = prep$sampInfo$event,
    k       = fixture$k,
    alpha   = 0.4,
    lambda  = 0.15,
    nu      = 0.5,
    lambdaW = 0.05,
    lambdaH = 0.05,
    verbose = FALSE,
    preprocess_info = preprocess_meta
  ))

  baseline <- predict(
    fit,
    newdata = desurv_data(prep$ex, prep$sampInfo$time, prep$sampInfo$event, fixture$k)
  )

  raw_subset <- fixture$X[prep$featInfo, , drop = FALSE]
  pred_raw <- predict(fit, newdata = raw_subset)

  target_values <- preprocess_meta$transform_target$values
  manual_norm <- preprocessCore::normalize.quantiles.use.target(
    raw_subset,
    target = target_values
  )
  if (!is.null(preprocess_meta$transform_target$genes)) {
    rownames(manual_norm) <- preprocess_meta$transform_target$genes
  }
  manual_pred <- predict(
    fit,
    newdata = desurv_data(manual_norm, prep$sampInfo$time, prep$sampInfo$event, fixture$k)
  )

  expect_equal(pred_raw, manual_pred)
  expect_equal(pred_raw, baseline)
})
