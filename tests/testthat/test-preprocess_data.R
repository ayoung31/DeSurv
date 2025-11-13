test_that("preprocess_data respects sample filters and supplied genes", {
  X <- matrix(
    1:12,
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:4))
  )
  y <- c(5, 4, 3, 2)
  d <- c(1, 0, 1, 0)
  dataset <- c("A", "A", "B", "B")

  res <- preprocess_data(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    samp_keeps         = c(1, 3, 4),
    genes              = c("g3", "g1"),
    method_trans_train = "none"
  )

  expect_equal(rownames(res$ex), c("g3", "g1"))
  expect_equal(colnames(res$ex), c("s1", "s3", "s4"))
  expect_equal(res$sampInfo$time, y[c(1, 3, 4)])
  expect_equal(res$samp_keeps, c(1, 3, 4))
})

test_that("preprocess_data fails when gene intersection is empty", {
  X <- matrix(
    c(10, 9, 1, 2,
      1, 2, 10, 9),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("g1", "g2"), paste0("s", 1:4))
  )
  y <- c(1, 2, 3, 4)
  d <- c(1, 1, 0, 0)
  dataset <- c("A", "A", "B", "B")

  expect_error(
    preprocess_data(
      X                  = X,
      y                  = y,
      d                  = d,
      dataset            = dataset,
      ngene              = 1,
      method_trans_train = "none"
    ),
    "No overlapping genes"
  )
})

test_that("desurv_cv preprocess option filters data and stores metadata", {
  X <- matrix(
    c(1, 2, 3,
      4, 5, 6),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )
  y <- c(5, 4, 3)
  d <- c(1, 0, 1)
  dataset <- c("A", "A", "B")
  samp_keeps <- c(TRUE, FALSE, TRUE)

  mock_cv_engine <- function(
      X, y, d, dataset,
      k_grid, alpha_grid, lambda_grid, nu_grid, lambdaW_grid, lambdaH_grid,
      ...) {
    expect_equal(dim(X), c(2, 2))
    expect_equal(y, c(5, 3))
    expect_equal(d, c(1, 1))
    expect_equal(dataset, c("A", "B"))

    cv_results <- data.frame(
      fold   = 1L,
      k      = k_grid[1],
      lambda = lambda_grid[1],
      nu     = nu_grid[1],
      lambdaW = lambdaW_grid[1],
      lambdaH = lambdaH_grid[1],
      init_id = 1L,
      alpha   = alpha_grid[1],
      cindex  = 0.5
    )

    summary_fold <- data.frame(
      fold = 1L,
      k = k_grid[1],
      lambda = lambda_grid[1],
      nu = nu_grid[1],
      lambdaW = lambdaW_grid[1],
      lambdaH = lambdaH_grid[1],
      alpha = alpha_grid[1],
      mean_cindex_fold = 0.5
    )

    summary <- data.frame(
      k = k_grid[1],
      lambda = lambda_grid[1],
      nu = nu_grid[1],
      lambdaW = lambdaW_grid[1],
      lambdaH = lambdaH_grid[1],
      alpha = alpha_grid[1],
      mean_cindex = 0.5,
      se_cindex = 0
    )

    list(
      cv_results = cv_results,
      summary_fold = summary_fold,
      summary = summary,
      alpha_grid = alpha_grid,
      hyper_grid = data.frame(
        k = k_grid,
        lambda = lambda_grid,
        nu = nu_grid,
        lambdaW = lambdaW_grid,
        lambdaH = lambdaH_grid
      ),
      folds = rep(1L, length(y))
    )
  }

  mock_fit <- function(...) {
    structure(list(cindex = 0.6), class = "desurv_fit")
  }

  testthat::local_mocked_bindings(
    .desurv_cv_grid_warmstart = mock_cv_engine,
    .desurv_cv_grid = mock_cv_engine,
    desurv_fit = mock_fit
  )

  fit <- desurv_cv(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    samp_keeps         = samp_keeps,
    genes              = rownames(X),
    method_trans_train = "none",
    preprocess         = TRUE,
    k_grid       = 2,
    alpha_grid   = 0.5,
    lambda_grid  = 1,
    nu_grid      = 0.5,
    lambdaW_grid = 0.1,
    lambdaH_grid = 0.1,
    n_starts     = 1,
    nfolds       = 2,
    tol          = 1e-4,
    maxit        = 10,
    verbose      = FALSE
  )

  expect_equal(fit$preprocess$samp_keeps, c(1, 3))
  expect_equal(fit$preprocess$genes, rownames(X))
})

test_that("desurv_cv_preprocess forwards arguments to desurv_cv", {
  called <- NULL
  mock_desurv_cv <- function(..., preprocess = FALSE) {
    called <<- list(...)
    expect_true(preprocess)
    structure(list(preprocess = list(tag = "ok")), class = "desurv_cv")
  }

  testthat::local_mocked_bindings(desurv_cv = mock_desurv_cv)

  expect_warning(
    res <- desurv_cv_preprocess(
      X = matrix(1:4, nrow = 2),
      y = c(1, 2),
      d = c(1, 0),
      dataset = c("A", "B"),
      genes = c("g1", "g2"),
      method_trans_train = "none",
      k_grid = 2,
      alpha_grid = 0.5,
      lambda_grid = 1,
      nu_grid = 0.5,
      lambdaW_grid = 0.1,
      lambdaH_grid = 0.1
    ),
    "Deprecated"
  )

  expect_s3_class(res, "desurv_cv")
  expect_equal(res$preprocess$tag, "ok")
})

test_that("preprocess_data verbose flag suppresses progress messages", {
  X <- matrix(
    runif(12),
    nrow = 3,
    dimnames = list(paste0("g", 1:3), paste0("s", 1:4))
  )
  y <- c(3, 2, 4, 5)
  d <- c(1, 1, 0, 1)
  dataset <- c("A", "A", "B", "B")

  expect_silent(
    preprocess_data(
      X                  = X,
      y                  = y,
      d                  = d,
      dataset            = dataset,
      ngene              = 2,
      method_trans_train = "none",
      verbose            = FALSE
    )
  )
})

test_that("preprocess_data rank transform produces column-wise ranks", {
  X <- matrix(
    c(5, 4, 3,
      2, 1, 0),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )
  y <- c(3, 2, 1)
  d <- c(1, 0, 1)
  dataset <- c("A", "A", "A")

  res <- preprocess_data(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    genes              = rownames(X),
    method_trans_train = "rank",
    verbose            = FALSE
  )

  expect_equal(unname(res$ex[1, ]), c(2, 2, 2))
  expect_equal(unname(res$ex[2, ]), c(1, 1, 1))
})

test_that("preprocess_data quant transform applies quantile normalization", {
  testthat::skip_if_not_installed("preprocessCore")
  X <- matrix(
    c(10, 2, 8,
      3,  7, 1),
    nrow = 2,
    dimnames = list(c("g1", "g2"), paste0("s", 1:3))
  )
  y <- c(2, 4, 6)
  d <- c(1, 0, 1)
  dataset <- rep("A", 3)

  res <- preprocess_data(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    genes              = rownames(X),
    method_trans_train = "quant",
    verbose            = FALSE
  )

  col1 <- res$ex[, 1]
  identical_cols <- vapply(
    seq_len(ncol(res$ex)),
    function(j) isTRUE(all.equal(res$ex[, j], col1)),
    logical(1)
  )
  expect_true(all(identical_cols))
})
