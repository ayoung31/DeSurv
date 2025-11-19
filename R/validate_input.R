# R/validate.R
.validate_desurv_data <- function(X, y, d, k, dataset = NULL) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be a numeric matrix.")
  if (anyNA(X) || !all(is.finite(X))) stop("X cannot contain NA/NaN/Inf.")

  p <- nrow(X); n <- ncol(X)
  if (p <= 1L || n <= 1L) stop("X must have at least 2 rows and 2 columns.")

  y <- as.numeric(y)
  d <- as.numeric(d)

  if (length(y) != n || length(d) != n) {
    stop("Lengths of y and d must equal ncol(X).")
  }
  if (anyNA(y) || anyNA(d) || !all(is.finite(y)) || !all(is.finite(d))) {
    stop("y and d cannot contain NA/NaN/Inf.")
  }
  if (!all(d %in% c(0,1))) stop("d must contain only 0s and 1s.")
  if (!all(y > 0)) stop("All survival times y must be > 0.")
  if (sum(d) <= 0) stop("There must be at least one event (sum(d) > 0).")

  k_num <- as.numeric(k)
  if (length(k_num) != 1L || is.na(k_num)) {
    stop("k must be a finite numeric scalar.")
  }
  k <- as.integer(k_num)
  if (k <= 0L || k > min(p, n)) {
    stop("k must be a positive integer <= min(nrow(X), ncol(X)).")
  }

  if (is.null(dataset)) {
    dataset <- rep("_all", n)
  } else {
    if (length(dataset) != n) {
      stop("`dataset` must have length equal to ncol(X).")
    }
    if (anyNA(dataset)) {
      stop("`dataset` cannot contain NA.")
    }
  }
  dataset <- as.character(dataset)

  list(
    X = X,
    y = y,
    d = d,
    k = k,
    dataset = dataset,
    p = p,
    n = n
  )
}

#' Construct a DeSurv dataset
#'
#' @description
#' Validate and encapsulate the core data required to fit a DeSurv model.
#' This constructor performs all structural checks on the inputs and returns
#' an object of class \code{"desurv_data"} that can be reused across multiple
#' model fits and cross-validation runs.
#'
#' @param X Numeric matrix (\eqn{p \times n}) of features (e.g., genes by samples).
#' @param y Numeric vector of survival times of length \eqn{n}.
#' @param d Numeric (0/1) vector of event indicators of length \eqn{n}.
#' @param k Positive integer rank (number of latent factors).
#' @param dataset Optional vector (length \code{ncol(X)}) identifying the
#'   dataset/cohort for each sample. Used for downstream stratified folds.
#'
#' @return
#' An object of class \code{"desurv_data"} containing the validated inputs.
#'
#' @export
desurv_data <- function(X, y, d, k, dataset = NULL) {
  args <- .validate_desurv_data(X, y, d, k, dataset = dataset)
  structure(args, class = "desurv_data")
}

#' @export
print.desurv_data <- function(x, ...) {
  cat("DeSurv dataset\n")
  cat("  p (features):", x$p, "\n")
  cat("  n (samples):", x$n, "\n")
  cat("  k (factors):", x$k, "\n")
  invisible(x)
}

# Internal: accept either raw inputs or a desurv_data object
.as_desurv_data <- function(X, y = NULL, d = NULL, k = NULL, dataset = NULL) {
  if (inherits(X, "desurv_data")) {
    if (!is.null(dataset)) {
      warning("Ignoring `dataset` because `X` is already a desurv_data object.")
    }
    return(X)
  }
  desurv_data(X = X, y = y, d = d, k = k, dataset = dataset)
}


.validate_desurv_hyperparams <- function(alpha, lambda, nu, lambdaW, lambdaH,
                                         tol, maxit, verbose) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha)) stop("alpha must be a scalar.")
  if (alpha < 0 || alpha >= 1) {
    warning("alpha outside [0,1); clipping.")
    alpha <- max(0, min(alpha, 1 - .Machine$double.eps))
  }

  nu <- as.numeric(nu)
  if (length(nu) != 1L || is.na(nu)) stop("nu must be a scalar.")
  if (nu < 0 || nu > 1) {
    warning("nu outside [0,1]; clipping.")
    nu <- max(0, min(nu, 1))
  }

  fix_pen <- function(v, name) {
    v <- as.numeric(v)
    if (length(v) != 1L || is.na(v)) stop(sprintf("%s must be scalar.", name))
    if (v < 0) {
      warning(sprintf("%s < 0; setting to 0.", name))
      v <- 0
    }
    v
  }

  lambda  <- fix_pen(lambda,  "lambda")
  lambdaW <- fix_pen(lambdaW, "lambdaW")
  lambdaH <- fix_pen(lambdaH, "lambdaH")

  tol <- as.numeric(tol)
  if (length(tol) != 1L || is.na(tol) || tol <= 0) {
    warning("tol not positive; setting tol = 1e-6.")
    tol <- 1e-6
  }

  maxit <- as.numeric(maxit)
  if (length(maxit) != 1L || is.na(maxit) || maxit <= 0) {
    warning("maxit not positive; setting maxit = 1000.")
    maxit <- 1000
  }
  maxit <- as.integer(maxit)

  verbose <- as.logical(verbose)
  if (length(verbose) != 1L || is.na(verbose)) {
    warning("verbose invalid; setting verbose = FALSE.")
    verbose <- FALSE
  }

  list(
    alpha   = alpha,
    lambda  = lambda,
    nu      = nu,
    lambdaW = lambdaW,
    lambdaH = lambdaH,
    tol     = tol,
    maxit   = maxit,
    verbose = verbose
  )
}

.validate_desurv_custom_init <- function(X, k, W0, H0, beta0) {
  if (is.null(W0)) {
    stop("`W0` must be supplied when providing custom initialization.", call. = FALSE)
  }

  W0 <- as.matrix(W0)
  if (!is.numeric(W0)) {
    stop("`W0` must be numeric.", call. = FALSE)
  }
  if (anyNA(W0) || !all(is.finite(W0))) {
    stop("`W0` cannot contain NA, NaN, or Inf.", call. = FALSE)
  }
  if (any(W0 < 0)) {
    stop("`W0` must be non-negative.", call. = FALSE)
  }

  p <- nrow(X)
  n <- ncol(X)
  if (nrow(W0) != p) stop("nrow(W0) must equal nrow(X).")
  if (ncol(W0) != k) stop("ncol(W0) must equal k.")

  if (is.null(H0)) {
    H0 <- .desurv_compute_nnls_H0(W0, X)
  } else {
    H0 <- as.matrix(H0)
    if (!is.numeric(H0)) stop("`H0` must be numeric.")
    if (anyNA(H0) || !all(is.finite(H0))) {
      stop("`H0` cannot contain NA, NaN, or Inf.")
    }
    if (any(H0 < 0)) stop("`H0` must be non-negative.")
  }

  if (nrow(H0) != k) stop("nrow(H0) must equal k.")
  if (ncol(H0) != n) stop("ncol(H0) must equal ncol(X).")

  if (is.null(beta0)) {
    beta0 <- rep(0, k)
  } else {
    beta0 <- as.numeric(beta0)
    if (length(beta0) != k) stop("length(beta0) must equal k.")
    if (anyNA(beta0) || !all(is.finite(beta0))) {
      stop("`beta0` cannot contain NA, NaN, or Inf.")
    }
  }

  list(W0 = W0, H0 = H0, beta0 = beta0)
}
