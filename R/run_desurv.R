#' Fit a DeSurv model from manuscript scripts
#'
#' Wrapper that preserves the legacy manuscript signature while delegating
#' to `DeSurv::desurv_fit()` under the hood.
run_desurv <- function(
    X, y, delta, k,
    alpha, lambda, eta,
    lambdaW = 0, lambdaH = 0,
    seed = NULL,
    tol = 1e-6,
    maxit = 1000,
    ninit = 20,
    imaxit = 30,
    W0 = NULL,
    H0 = NULL,
    beta0 = NULL,
    verbose = TRUE
) {
  if (!requireNamespace("DeSurv", quietly = TRUE)) {
    stop("Package `DeSurv` must be installed.")
  }

  data_obj <- DeSurv::desurv_data(X, y, delta, k)
  tol_init <- max(tol * 10, 1e-4)

  DeSurv::desurv_fit(
    data_obj,
    alpha         = alpha,
    lambda        = lambda,
    nu            = eta,
    lambdaW       = lambdaW,
    lambdaH       = lambdaH,
    W0            = W0,
    H0            = H0,
    beta0         = beta0,
    seed          = seed,
    tol           = tol,
    tol_init      = tol_init,
    maxit         = maxit,
    imaxit        = imaxit,
    ninit         = ninit,
    parallel_init = FALSE,
    verbose       = verbose
  )
}
