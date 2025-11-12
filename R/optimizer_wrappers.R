#' @keywords internal
.run_optimize_loss <- function(
    X,
    y,
    d,
    W0,
    H0,
    beta0,
    alpha,
    lambda,
    nu,
    lambdaW,
    lambdaH,
    tol,
    maxit,
    verbose,
    theta_init        = 0.5,
    rho               = 0.5,
    max_backtracks    = 10L,
    eps_beta          = 1e-6,
    max_iter_beta     = 100L,
    max_iter_beta_final = 500L
) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  d <- as.numeric(d)

  storage.mode(X) <- "double"
  storage.mode(y) <- "double"
  storage.mode(d) <- "double"

  k <- ncol(W0)
  n <- ncol(X)
  p <- nrow(X)
  n_event <- sum(d)

  if (n_event <= 0) {
    stop("At least one event (d == 1) is required to fit the survival model.")
  }

  theta_init <- as.numeric(theta_init)
  rho <- as.numeric(rho)
  max_backtracks <- as.integer(max_backtracks)
  eps_beta <- as.numeric(eps_beta)
  max_iter_beta <- as.integer(max_iter_beta)
  max_iter_beta_final <- as.integer(max_iter_beta_final)

  # C++ optimizer should stay silent; verbose controls only R-side messaging.
  cpp_verbose <- FALSE

  optimize_loss_cpp(
    X, y, d,
    k, n, p, n_event,
    W0, H0, beta0,
    alpha, lambda, nu,
    lambdaW, lambdaH,
    tol, maxit, cpp_verbose,
    theta_init, rho, max_backtracks,
    eps_beta, max_iter_beta, max_iter_beta_final
  )
}
