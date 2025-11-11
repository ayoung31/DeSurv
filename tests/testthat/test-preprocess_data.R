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

test_that("desurv_cv_preprocess forwards processed data and metadata", {
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

  mock_cv <- function(X, y, d, ...) {
    expect_equal(dim(X), c(2, 2))
    expect_equal(y, c(5, 3))
    expect_equal(d, c(1, 1))
    list(fit = "mock")
  }

  testthat::local_mocked_bindings(desurv_cv = mock_cv)

  fit <- desurv_cv_preprocess(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    samp_keeps         = samp_keeps,
    genes              = rownames(X),
    method_trans_train = "none"
  )

  expect_equal(fit$preprocess$samp_keeps, c(1, 3))
})
