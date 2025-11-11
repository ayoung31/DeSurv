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
