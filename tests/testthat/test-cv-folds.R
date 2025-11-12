test_that("folds stratify by event and dataset", {
  set.seed(321)
  n <- 30
  d <- sample(c(0, 1), size = n, replace = TRUE)
  dataset <- rep(c("A", "B", "C"), length.out = n)

  folds <- DeSurv:::`.desurv_make_folds_stratified`(
    d,
    dataset = dataset,
    nfolds  = 5,
    seed    = 99
  )

  expect_length(folds, n)
  strata <- interaction(d, dataset, drop = TRUE)
  tab <- table(strata, folds)
  balance_ok <- apply(tab, 1, function(x) max(x) - min(x) <= 1)
  expect_true(all(balance_ok))
})
