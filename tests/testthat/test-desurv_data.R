test_that("desurv_data validates core inputs", {
  fixture <- make_fixture_dataset()
  dd <- desurv_data(fixture$X, fixture$y, fixture$d, fixture$k)

  expect_s3_class(dd, "desurv_data")
  expect_equal(dd$k, fixture$k)
  expect_equal(dd$p, fixture$p)
  expect_equal(dd$n, fixture$n)

  bad_y <- fixture$y[-1]
  expect_error(
    desurv_data(fixture$X, bad_y, fixture$d, fixture$k),
    "Lengths of y and d must equal ncol",
    fixed = FALSE
  )

  bad_d <- fixture$d
  bad_d[] <- 0
  expect_error(
    desurv_data(fixture$X, fixture$y, bad_d, fixture$k),
    "at least one event",
    fixed = FALSE
  )
})

test_that(".validate_desurv_custom_init enforces dimensions and positivity", {
  fixture <- make_fixture_dataset()
  init <- .validate_desurv_custom_init(
    fixture$X,
    fixture$k,
    W0 = matrix(runif(fixture$p * fixture$k), nrow = fixture$p),
    H0 = matrix(runif(fixture$k * fixture$n), nrow = fixture$k),
    beta0 = rep(0, fixture$k)
  )
  expect_equal(dim(init$W0), c(fixture$p, fixture$k))

  H_bad <- matrix(runif(fixture$k * fixture$n), nrow = fixture$k + 1)
  expect_error(
    .validate_desurv_custom_init(
      fixture$X,
      fixture$k,
      W0 = init$W0,
      H0 = H_bad,
      beta0 = init$beta0
    ),
    "nrow\\(H0\\) must equal k"
  )

  expect_error(
    .validate_desurv_custom_init(
      fixture$X,
      fixture$k,
      W0 = init$W0 - 10,
      H0 = init$H0,
      beta0 = init$beta0
    ),
    "must be non-negative"
  )
})


test_that("desurv_data stores dataset labels", {
  fixture <- make_fixture_dataset()
  dataset <- rep(c("A", "B"), length.out = fixture$n)
  dd <- desurv_data(
    fixture$X,
    fixture$y,
    fixture$d,
    fixture$k,
    dataset = dataset
  )
  expect_equal(dd$dataset, as.character(dataset))
})

# Tests for Issue #1: Missing Nonnegativity Validation for X
test_that("desurv_data rejects negative X values (NMF requirement)", {
  fixture <- make_fixture_dataset()

  # Create X with some negative values
  X_negative <- fixture$X
  X_negative[1, 1] <- -1

  expect_error(
    desurv_data(X_negative, fixture$y, fixture$d, fixture$k),
    "non-negative",
    fixed = FALSE
  )

  # All negative should also fail
  X_all_negative <- -abs(fixture$X)
  expect_error(
    desurv_data(X_all_negative, fixture$y, fixture$d, fixture$k),
    "non-negative",
    fixed = FALSE
  )
})

# Tests for Issue #5: All-Zero X Yields Division by Zero
test_that("desurv_data rejects all-zero X (Xnorm=0)",
{
  fixture <- make_fixture_dataset()

  # All zeros
  X_zeros <- matrix(0, nrow = fixture$p, ncol = fixture$n)
  expect_error(
    desurv_data(X_zeros, fixture$y, fixture$d, fixture$k),
    "cannot be all zeros",
    fixed = FALSE
  )
})

# Test that X with at least one non-zero value is accepted
test_that("desurv_data accepts X with at least one non-zero value", {
  fixture <- make_fixture_dataset()

  # Nearly all zeros but one non-zero value
  X_sparse <- matrix(0, nrow = fixture$p, ncol = fixture$n)
  X_sparse[1, 1] <- 1

  dd <- desurv_data(X_sparse, fixture$y, fixture$d, fixture$k)
  expect_s3_class(dd, "desurv_data")
})

# Tests for Issue #6: Alpha Clipping Contradicts Documentation
test_that("alpha=1 is clipped with warning", {
  # Note: alpha is validated in .validate_desurv_hyperparams, not desurv_data
  # Test the hyperparameter validation directly
  expect_warning(
    result <- .validate_desurv_hyperparams(
      alpha = 1,
      lambda = 0.1,
      nu = 0.5,
      lambdaW = 0,
      lambdaH = 0,
      tol = 1e-6,
      maxit = 100,
      verbose = FALSE
    ),
    "alpha outside"
  )
  expect_true(result$alpha < 1)
  expect_true(result$alpha > 0.99)  # Should be very close to 1

  # alpha > 1 should also be clipped
  expect_warning(
    result2 <- .validate_desurv_hyperparams(
      alpha = 1.5,
      lambda = 0.1,
      nu = 0.5,
      lambdaW = 0,
      lambdaH = 0,
      tol = 1e-6,
      maxit = 100,
      verbose = FALSE
    ),
    "alpha outside"
  )
  expect_true(result2$alpha < 1)
})
