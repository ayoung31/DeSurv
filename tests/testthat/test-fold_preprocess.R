fold_helper <- getFromNamespace(".desurv_prepare_fold_payload", "DeSurv")

test_that("fold preprocessing filters genes per fold", {
  set.seed(123)
  X <- matrix(runif(5 * 8), nrow = 5, ncol = 8)
  rownames(X) <- paste0("g", seq_len(5))
  y <- rexp(8)
  d <- sample(0:1, 8, replace = TRUE)
  dataset <- rep(c("A", "B"), each = 4)

  payload <- fold_helper(
    X_full = X,
    y_full = y,
    d_full = d,
    dataset_full = dataset,
    idx_tr = 1:6,
    idx_val = 7:8,
    preprocess = TRUE,
    ngene = 3,
    genes = NULL,
    method_trans_train = "rank",
    verbose = FALSE
  )

  expect_true(nrow(payload$X_tr) <= 3)
  expect_equal(rownames(payload$X_tr), rownames(payload$X_val))
  expect_true(all(apply(payload$X_val, 2, function(col) all(col %in% seq_len(nrow(payload$X_tr))))))
})

test_that("fold preprocessing handles raw inputs when disabled", {
  set.seed(321)
  X <- matrix(runif(4 * 6), nrow = 4, ncol = 6)
  rownames(X) <- paste0("g", seq_len(4))
  y <- rexp(6)
  d <- sample(0:1, 6, replace = TRUE)
  dataset <- rep("A", 6)

  payload <- fold_helper(
    X_full = X,
    y_full = y,
    d_full = d,
    dataset_full = dataset,
    idx_tr = 1:4,
    idx_val = 5:6,
    preprocess = FALSE,
    ngene = 2,
    genes = NULL,
    method_trans_train = "rank",
    verbose = FALSE
  )

  expect_equal(nrow(payload$X_tr), nrow(X))
  expect_equal(rownames(payload$X_tr), rownames(X))
  expect_equal(rownames(payload$X_val), rownames(X))
})
