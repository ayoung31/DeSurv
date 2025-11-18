test_that("desurv_get_top_genes returns structured outputs", {
  W <- matrix(
    c(
      5, 1,
      4, 0.5,
      0.2, 3,
      0.1, 2
    ),
    nrow = 4,
    ncol = 2,
    byrow = TRUE
  )
  rownames(W) <- paste0("g", seq_len(4))
  colnames(W) <- c("f1", "f2")

  res <- desurv_get_top_genes(W, ntop = 2)
  expect_type(res, "list")
  expect_s3_class(res$top_genes, "data.frame")
  expect_equal(nrow(res$top_genes), 2)
  expect_equal(colnames(res$top_genes), colnames(W))
  expect_equal(res$top_genes[1, "f1"], "g1")
  expect_equal(res$top_genes[1, "f2"], "g3")
  expect_s3_class(res$top_indices, "data.frame")
  expect_equal(res$top_indices[1, "f1"], 1)
  expect_s3_class(res$diffs, "data.frame")
})

test_that("desurv_cv records ntop in best configuration", {
  skip_on_cran()
  fixture <- make_fixture_dataset(n = 12, seed = 42)

  fit_ntop <- desurv_cv(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    k_grid = 2,
    alpha_grid = c(0.2, 0.4),
    lambda_grid = 0.1,
    nu_grid = 0.5,
    lambdaW_grid = 0,
    lambdaH_grid = 0,
    n_starts = 1,
    nfolds = 2,
    ntop = 2,
    seed = 99,
    engine = "cold",
    verbose = FALSE
  )

  expect_equal(fit_ntop$best$ntop, 2)
  expect_equal(fit_ntop$ntop, 2)

  fit_no_ntop <- desurv_cv(
    X = fixture$X,
    y = fixture$y,
    d = fixture$d,
    k_grid = 2,
    alpha_grid = c(0.2, 0.4),
    lambda_grid = 0.1,
    nu_grid = 0.5,
    lambdaW_grid = 0,
    lambdaH_grid = 0,
    n_starts = 1,
    nfolds = 2,
    ntop = NULL,
    seed = 100,
    engine = "cold",
    verbose = FALSE
  )

  expect_null(fit_no_ntop$best$ntop)
  expect_null(fit_no_ntop$ntop)
})
