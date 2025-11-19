test_that("desurv_consensus_seed builds consensus W0/H0", {
  fixture <- make_fixture_dataset(p = 6, n = 10, k = 2)

  make_fit <- function(scale_vec) {
    W <- cbind(scale_vec, rev(scale_vec))
    rownames(W) <- paste0("gene", seq_len(fixture$p))
    structure(
      list(
        W = W,
        H = matrix(runif(fixture$k * fixture$n), nrow = fixture$k),
        beta = runif(fixture$k),
        cindex = runif(1)
      ),
      class = "desurv_fit"
    )
  }

  fits <- list(
    make_fit(c(5, 4, 3, 2, 1, 0.5)),
    make_fit(c(4, 5, 2, 1, 3, 0.2)),
    make_fit(c(3, 1, 5, 4, 2, 0.1))
  )

  consensus <- desurv_consensus_seed(
    fits = fits,
    X = fixture$X,
    ntop = 3,
    clustering = "average"
  )

  expect_true(is.matrix(consensus$W0))
  expect_true(is.matrix(consensus$H0))
  expect_equal(dim(consensus$W0), c(fixture$p, fixture$k))
  expect_equal(dim(consensus$H0), c(fixture$k, fixture$n))
  expect_true(all(consensus$W0 >= 0))
  expect_true(all(consensus$H0 >= 0))
  expect_true(is.matrix(consensus$consensus))
  expect_equal(dim(consensus$consensus), c(fixture$p, fixture$p))
  expect_true(consensus$n_factors > 0)
})

test_that("custom initialization infers H0 via NNLS when omitted", {
  fixture <- make_fixture_dataset(p = 5, n = 6, k = 2)
  W0 <- matrix(runif(fixture$p * fixture$k, min = 0, max = 2),
               nrow = fixture$p, ncol = fixture$k)

  expect_silent({
    init_vals <- DeSurv:::.validate_desurv_custom_init(
      fixture$X,
      fixture$k,
      W0 = W0,
      H0 = NULL,
      beta0 = NULL
    )
    expect_equal(dim(init_vals$H0), c(fixture$k, fixture$n))
    expect_true(all(init_vals$H0 >= 0))
    expect_equal(length(init_vals$beta0), fixture$k)
  })
})
