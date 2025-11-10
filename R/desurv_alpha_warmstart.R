#' Warm-start DeSurv over alpha for multiple random initializations (optional parallel)
#'
#' @description
#' Run DeSurv fits over a sequence of supervision parameters \code{alpha_grid}
#' using warm starts, and repeat this for multiple random initializations.
#' For each initialization path, a single random initialization
#' \eqn{(W_0, H_0, \beta_0)} is drawn, and a full warmstart trajectory is fit
#' over \code{alpha_grid}. Paths can optionally be executed in parallel.
#'
#' @inheritParams desurv_fit
#' @param alpha_grid Numeric vector of supervision parameter values to explore.
#' @param n_starts Integer number of independent random initialization paths.
#' @param seed Integer or \code{NULL}; base seed for reproducibility.
#' @param tol,maxit Numeric; optimizer tolerance and max iterations per fit.
#' @param parallel Logical; if \code{TRUE}, run different initialization paths
#'   in parallel using \code{parallel::mclapply()} (Unix/macOS only).
#' @param ncores Integer or \code{NULL}; number of cores to use when
#'   \code{parallel = TRUE}. Defaults to \code{parallel::detectCores()} capped
#'   at \code{n_starts}.
#' @param verbose Logical; print progress.
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{results}: data frame with \code{init_id}, \code{alpha}, \code{cindex}.
#'   \item \code{fits}: list of all \code{"desurv_fit"} objects by path and alpha.
#'   \item \code{alpha_grid}: vector of alpha values used.
#' }
#'
#' @export
desurv_alpha_warmstart <- function(
    X, y = NULL, d = NULL, k = NULL,
    alpha_grid,
    lambda, nu, lambdaW, lambdaH,
    n_starts = 5,
    seed     = 123,
    tol      = 1e-6,
    maxit    = 1000,
    parallel = FALSE,
    ncores   = NULL,
    verbose  = TRUE
) {
  # Validate alpha grid
  alpha_grid <- as.numeric(alpha_grid)
  if (length(alpha_grid) == 0L || anyNA(alpha_grid))
    stop("alpha_grid must be a non-empty numeric vector without NA.")
  A <- length(alpha_grid)

  # Validate n_starts
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts <= 0)
    stop("n_starts must be a positive integer.")

  # Validate and prepare data
  data_obj <- .as_desurv_data(X, y = y, d = d, k = k)
  p <- data_obj$p; n <- data_obj$n; k <- data_obj$k
  max_x <- max(data_obj$X)

  # Internal worker function: run one initialization path
  run_one_path <- function(j) {
    if (!is.null(seed)) set.seed(as.integer(seed + j - 1L))
    if (verbose && !parallel)
      message(sprintf("Warmstart path %d / %d", j, n_starts))

    # Draw random initialization
    W0 <- matrix(runif(p * k, 0, max_x), nrow = p, ncol = k)
    H0 <- matrix(runif(n * k, 0, max_x), nrow = k, ncol = n)
    beta0 <- rep(0, k)

    fits_j <- vector("list", A)
    cindex_j <- numeric(A)

    for (i in seq_len(A)) {
      a_i <- alpha_grid[i]

      fit_i <- desurv_fit(
        X       = data_obj,
        alpha   = a_i,
        lambda  = lambda,
        nu      = nu,
        lambdaW = lambdaW,
        lambdaH = lambdaH,
        W0      = W0,
        H0      = H0,
        beta0   = beta0,
        tol     = tol,
        maxit   = maxit,
        verbose = verbose && !parallel
      )

      fits_j[[i]] <- fit_i
      cindex_j[i] <- fit_i$cindex

      # Warmstart for next alpha
      W0 <- fit_i$W; H0 <- fit_i$H; beta0 <- fit_i$beta
    }

    list(fits = fits_j, cindex = cindex_j)
  }

  # Run either sequentially or in parallel
  if (parallel) {
    if (.Platform$OS.type == "windows") {
      warning("Parallel warmstart is not supported on Windows; falling back to sequential.")
      parallel <- FALSE
    }
  }

  if (parallel) {
    if (is.null(ncores)) ncores <- parallel::detectCores()
    ncores <- min(ncores, n_starts)
    if (verbose)
      message(sprintf("Running %d warmstart paths in parallel on %d cores.", n_starts, ncores))

    paths <- parallel::mclapply(
      seq_len(n_starts),
      run_one_path,
      mc.cores = ncores
    )
  } else {
    paths <- lapply(seq_len(n_starts), run_one_path)
  }

  # Aggregate results
  fits <- lapply(paths, `[[`, "fits")
  cmat <- do.call(rbind, lapply(paths, `[[`, "cindex"))

  results <- data.frame(
    init_id = rep(seq_len(n_starts), each = A),
    alpha   = rep(alpha_grid, times = n_starts),
    cindex  = as.vector(t(cmat)),
    row.names = NULL
  )

  list(
    results    = results,
    fits       = fits,
    alpha_grid = alpha_grid
  )
}
