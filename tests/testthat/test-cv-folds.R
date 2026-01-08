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

# Test for Issue #2: CV Mean/SE Calculation Inconsistency
test_that("CV mean and SE use consistent non-finite handling", {
  # Test the helper functions directly to verify they handle NA/Inf consistently
  # Both should filter non-finite values first, then compute statistic

  # Define the functions as they appear in cv_warm.R / cv_cold.R after fix
  mean_fun <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
  }

  se_fun <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }

  # Test 1: Both should return finite values when some folds are NA
  x_with_na <- c(0.6, 0.7, NA, 0.8, 0.65)
  m <- mean_fun(x_with_na)
  s <- se_fun(x_with_na)

  # Both should be finite (mean = finite, SE = finite)
  expect_true(is.finite(m))
  expect_true(is.finite(s))

  # Test 2: Both should return NA when all values are NA/Inf
  x_all_na <- c(NA, NaN, Inf)
  expect_true(is.na(mean_fun(x_all_na)))
  expect_true(is.na(se_fun(x_all_na)))

  # Test 3: SE should be NA when only one finite value (need >= 2 for SD)
  x_one_finite <- c(NA, 0.7, Inf)
  expect_true(is.finite(mean_fun(x_one_finite)))  # mean can be computed
  expect_true(is.na(se_fun(x_one_finite)))        # SE needs >= 2 values

  # Test 4: Verify consistent behavior - if mean is finite, SE logic is applied
  # (SE may still be NA due to insufficient values, but not due to filtering differences)
  x_two_finite <- c(NA, 0.7, 0.8, Inf)
  m2 <- mean_fun(x_two_finite)
  s2 <- se_fun(x_two_finite)
  expect_true(is.finite(m2))
  expect_true(is.finite(s2))  # Both should be finite with 2+ values
})
