#' @keywords internal
#' @title Multi-start initialization for DeSurv
#' @description
#' Run \code{ninit} short DeSurv optimizations with at most \code{imaxit}
#' iterations each, using a relaxed tolerance \code{tol_init}, evaluate
#' C-index of each, and return the best initialization (W, H, beta).
#' Can optionally run the short initializations in parallel.
init <- function(
    X, y, delta, k,
    alpha, lambda, nu, lambdaW, lambdaH,
    seed     = 123,
    tol_init = 1e-4,
    imaxit   = 30,
    verbose  = TRUE,
    ninit    = 20,
    parallel = FALSE,
    ncores   = NULL
) {
  p <- nrow(X)
  n <- ncol(X)

  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  # helper to run ONE short initialization
  run_one_init <- function(i, verbose_inner = FALSE) {
    if (verbose_inner) {
      message(sprintf("  Short initialization %d / %d", i, ninit))
    }

    max_x <- max(X)
    W0_i    <- matrix(runif(p * k, 0, max_x), nrow = p, ncol = k)
    H0_i    <- matrix(runif(n * k, 0, max_x), nrow = k, ncol = n)
    beta0_i <- rep(0, k)

    fit_short <- tryCatch(
      optimize_loss_cpp(
        X, y, delta,
        W0_i, H0_i, beta0_i,
        alpha, lambda, nu,
        lambdaW, lambdaH,
        tol_init, imaxit, verbose_inner, TRUE
      ),
      error = function(e) {
        if (verbose_inner) {
          warning(sprintf("  Initialization %d failed in optimize_loss_cpp(): %s",
                          i, e$message))
        }
        return(NULL)
      }
    )

    if (is.null(fit_short) || isTRUE(fit_short$nan_flag)) {
      return(NULL)
    }

    lp_short <- tryCatch({
      Z_short <- t(X) %*% fit_short$W
      drop(Z_short %*% fit_short$beta)
    }, error = function(e) {
      if (verbose_inner) {
        warning(sprintf("  Initialization %d failed computing lp: %s",
                        i, e$message))
      }
      return(rep(NA_real_, n))
    })

    if (anyNA(lp_short) || !all(is.finite(lp_short))) {
      return(NULL)
    }

    cval_short <- tryCatch(
      cvwrapr::getCindex(lp_short, survival::Surv(y, delta)),
      error = function(e) {
        if (verbose_inner) {
          warning(sprintf("  Initialization %d failed in getCindex(): %s",
                          i, e$message))
        }
        return(NA_real_)
      }
    )

    if (is.na(cval_short)) {
      return(NULL)
    }

    list(
      cindex = cval_short,
      W      = fit_short$W,
      H      = fit_short$H,
      beta   = fit_short$beta
    )
  }

  if (verbose) {
    message(sprintf(
      "Running %d short initializations (imaxit = %d, tol_init = %g)%s.",
      ninit, imaxit, tol_init,
      if (parallel) " in parallel" else ""
    ))
  }

  results <- vector("list", ninit)

  if (parallel) {
    # fall back to sequential on Windows
    if (.Platform$OS.type == "windows") {
      if (verbose) {
        warning("Parallel initialization is not supported on Windows; falling back to sequential.")
      }
      for (i in seq_len(ninit)) {
        results[[i]] <- run_one_init(i, verbose_inner = verbose)
      }
    } else {
      if (is.null(ncores)) {
        ncores <- parallel::detectCores()
      }
      ncores <- max(1L, min(as.integer(ncores), ninit))

      # In parallel, suppress per-init verbose to avoid messy interleaving
      results <- parallel::mclapply(
        seq_len(ninit),
        function(i) run_one_init(i, verbose_inner = FALSE),
        mc.cores = ncores
      )
    }
  } else {
    # sequential
    for (i in seq_len(ninit)) {
      results[[i]] <- run_one_init(i, verbose_inner = verbose)
    }
  }

  # pick best non-NULL result
  cvals <- sapply(results, function(res) if (is.null(res)) NA_real_ else res$cindex)
  if (all(is.na(cvals))) {
    stop("All short initializations failed or produced invalid fits.")
  }

  best_idx <- which.max(cvals)
  best_res <- results[[best_idx]]

  if (verbose) {
    message(sprintf("Best short initialization C-index: %.4f (init %d).",
                    best_res$cindex, best_idx))
  }

  list(
    W0          = best_res$W,
    H0          = best_res$H,
    beta0       = best_res$beta,
    cindex_init = best_res$cindex
  )
}
