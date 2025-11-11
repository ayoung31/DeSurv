#' @keywords internal
#' @title Stratified fold assignment by event indicator
#' @description
#' Create stratified K-fold assignments so that events (d == 1) and
#' non-events (d == 0) are approximately balanced across folds.
.desurv_make_folds_stratified <- function(d, nfolds, seed = NULL) {
  d <- as.numeric(d)
  if (any(is.na(d)) || !all(d %in% c(0, 1))) {
    stop("d must be a 0/1 vector without NA for stratified folds.")
  }

  nfolds <- as.integer(nfolds)
  if (length(nfolds) != 1L || is.na(nfolds) || nfolds <= 1L) {
    stop("nfolds must be an integer >= 2.")
  }

  n <- length(d)
  if (!is.null(seed)) set.seed(as.integer(seed))

  idx_event <- which(d == 1)
  idx_cens  <- which(d == 0)

  # shuffle within strata
  idx_event <- sample(idx_event, length(idx_event))
  idx_cens  <- sample(idx_cens,  length(idx_cens))

  # assign folds within each stratum
  folds_event <- rep(seq_len(nfolds), length.out = length(idx_event))
  folds_cens  <- rep(seq_len(nfolds), length.out = length(idx_cens))

  folds <- integer(n)
  folds[idx_event] <- folds_event
  folds[idx_cens]  <- folds_cens

  # optional: warn if some fold has no events
  events_per_fold <- tapply(d, folds, sum)
  if (any(events_per_fold == 0)) {
    warning("Some folds have zero events; this is unavoidable when events < nfolds.")
  }

  folds
}



#' @keywords internal
#' @title Build and validate DeSurv hyperparameter grid
.desurv_make_hyper_grid <- function(
    k_grid,
    lambda_grid,
    nu_grid,
    lambdaW_grid,
    lambdaH_grid
) {
  k_grid       <- unique(as.integer(k_grid))
  lambda_grid  <- as.numeric(lambda_grid)
  nu_grid      <- as.numeric(nu_grid)
  lambdaW_grid <- as.numeric(lambdaW_grid)
  lambdaH_grid <- as.numeric(lambdaH_grid)

  if (length(k_grid) == 0L)       stop("k_grid must have length > 0.")
  if (length(lambda_grid) == 0L)  stop("lambda_grid must have length > 0.")
  if (length(nu_grid) == 0L)      stop("nu_grid must have length > 0.")
  if (length(lambdaW_grid) == 0L) stop("lambdaW_grid must have length > 0.")
  if (length(lambdaH_grid) == 0L) stop("lambdaH_grid must have length > 0.")

  hyper_grid <- expand.grid(
    k       = k_grid,
    lambda  = lambda_grid,
    nu      = nu_grid,
    lambdaW = lambdaW_grid,
    lambdaH = lambdaH_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  list(
    hyper_grid  = hyper_grid,
    k_grid      = k_grid,
    lambda_grid = lambda_grid,
    nu_grid     = nu_grid,
    lambdaW_grid = lambdaW_grid,
    lambdaH_grid = lambdaH_grid
  )
}

#' @keywords internal
#' @title CV over hyperparameter grid with alpha warmstarts
#' @description
#' Internal helper for cross-validation over a grid of hyperparameters
#' \code{k, alpha, lambda, nu, lambdaW, lambdaH}, using warmstarts along
#' \code{alpha_grid} via \code{desurv_alpha_warmstart()} for each combination
#' of \code{k, lambda, nu, lambdaW, lambdaH}. Parallelism (if enabled) is
#' applied at the (hyperparameter, fold) job level; within each job,
#' \code{desurv_alpha_warmstart} is always called with \code{parallel = FALSE}
#' to avoid nested parallelism.
#' @export
.desurv_cv_grid_warmstart <- function(
    X, y = NULL, d = NULL,
    k_grid,
    alpha_grid,
    lambda_grid,
    nu_grid,
    lambdaW_grid,
    lambdaH_grid,
    n_starts      = 3,
    nfolds        = 5,
    folds         = NULL,
    seed          = 123,
    tol           = 1e-6,
    maxit         = 1000,
    parallel_grid = FALSE,
    ncores_grid   = NULL,
    verbose       = TRUE
) {
  ## --- alpha grid ---
  alpha_grid <- as.numeric(alpha_grid)
  if (length(alpha_grid) == 0L || anyNA(alpha_grid)) {
    stop("alpha_grid must be a non-empty numeric vector without NA.")
  }
  A <- length(alpha_grid)

  ## --- hyperparameter grid ---
  hg <- .desurv_make_hyper_grid(
    k_grid       = k_grid,
    lambda_grid  = lambda_grid,
    nu_grid      = nu_grid,
    lambdaW_grid = lambdaW_grid,
    lambdaH_grid = lambdaH_grid
  )
  hyper_grid <- hg$hyper_grid
  H <- nrow(hyper_grid)

  ## --- n_starts / nfolds checks ---
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts <= 0L) {
    stop("n_starts must be a positive integer.")
  }

  nfolds <- as.integer(nfolds)
  if (length(nfolds) != 1L || is.na(nfolds) || nfolds <= 1L) {
    stop("nfolds must be an integer >= 2.")
  }

  ## --- validate data once ---
  if (inherits(X, "desurv_data")) {
    X_full <- X$X
    y_full <- X$y
    d_full <- X$d
  } else {
    X_full <- as.matrix(X)
    y_full <- y
    d_full <- d
    tmp    <- .validate_desurv_data(X_full, y_full, d_full, k = hyper_grid$k[1L])
    X_full <- tmp$X; y_full <- tmp$y; d_full <- tmp$d
  }

  n <- ncol(X_full)

  ## --- folds: stratified by event ---
  if (is.null(folds)) {
    folds <- .desurv_make_folds_stratified(d_full, nfolds = nfolds, seed = seed)
  } else {
    folds <- as.integer(folds)
    if (length(folds) != n) {
      stop("folds must have length equal to ncol(X).")
    }
    if (!all(folds %in% seq_len(nfolds))) {
      stop("All entries in folds must be in 1:nfolds.")
    }
  }

  if (verbose) {
    message(sprintf(
      "Running CV over %d hyperparameter combos, %d alphas, %d starts, %d folds.",
      H, A, n_starts, nfolds
    ))
  }

  ## --- define (hyper, fold) jobs ---
  jobs <- expand.grid(
    h = seq_len(H),
    f = seq_len(nfolds),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  ## --- inner job function ---
  run_one_job <- function(job_row) {
    h <- job_row$h
    f <- job_row$f

    hp <- hyper_grid[h, ]
    k_cur       <- hp$k
    lambda_cur  <- hp$lambda
    nu_cur      <- hp$nu
    lambdaW_cur <- hp$lambdaW
    lambdaH_cur <- hp$lambdaH

    if (verbose && !parallel_grid) {
      message(sprintf(
        "Hyper combo %d/%d (k=%d, lambda=%.3g, nu=%.3g, lambdaW=%.3g, lambdaH=%.3g), fold %d/%d",
        h, H, k_cur, lambda_cur, nu_cur, lambdaW_cur, lambdaH_cur, f, nfolds
      ))
    }

    idx_val <- which(folds == f)
    idx_tr  <- setdiff(seq_len(n), idx_val)

    X_tr <- X_full[, idx_tr, drop = FALSE]
    y_tr <- y_full[idx_tr]
    d_tr <- d_full[idx_tr]

    X_val <- X_full[, idx_val, drop = FALSE]
    y_val <- y_full[idx_val]
    d_val <- d_full[idx_val]

    data_tr <- desurv_data(X_tr, y_tr, d_tr, k_cur)

    # unique seed for this (hyper, fold)
    seed_fold <- if (is.null(seed)) NULL else as.integer(seed + 10000L * (h - 1L) + 100L * (f - 1L))

    # warmstart paths on TRAINING data for this combo (sequential inside)
    ws_res <- desurv_alpha_warmstart(
      X          = data_tr,
      alpha_grid = alpha_grid,
      lambda     = lambda_cur,
      nu         = nu_cur,
      lambdaW    = lambdaW_cur,
      lambdaH    = lambdaH_cur,
      n_starts   = n_starts,
      seed       = seed_fold,
      tol        = tol,
      maxit      = maxit,
      parallel   = FALSE,    # important: avoid nested parallelism
      ncores     = NULL,
      verbose    = verbose && !parallel_grid
    )

    # compute validation C-index matrix: n_starts x A
    cmat_val <- matrix(NA_real_, nrow = n_starts, ncol = A)
    for (init_id in seq_len(n_starts)) {
      for (i in seq_len(A)) {
        fit_ji <- ws_res$fits[[init_id]][[i]]
        lp_val <- predict(fit_ji, newdata = X_val, type = "lp")

        cmat_val[init_id, i] <-
          if (sum(d_val) == 0L) NA_real_
        else cvwrapr::getCindex(lp_val, survival::Surv(y_val, d_val))
      }
    }

    # tidy data.frame for this (hyper, fold) job
    df <- expand.grid(
      init_id = seq_len(n_starts),
      alpha   = alpha_grid,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
    df$cindex <- as.vector(t(cmat_val))
    df$fold   <- f
    df$k       <- k_cur
    df$lambda  <- lambda_cur
    df$nu      <- nu_cur
    df$lambdaW <- lambdaW_cur
    df$lambdaH <- lambdaH_cur

    # reorder columns a bit
    df <- df[, c("fold", "k", "lambda", "nu", "lambdaW", "lambdaH",
                 "init_id", "alpha", "cindex")]
    df
  }


  ## --- run jobs (optionally in parallel) ---
  if (parallel_grid && .Platform$OS.type == "windows") {
    warning("parallel_grid = TRUE is not supported on Windows; falling back to sequential.")
    parallel_grid <- FALSE
  }

  if (parallel_grid) {
    if (is.null(ncores_grid)) {
      ncores_grid <- parallel::detectCores()
    }
    ncores_grid <- max(1L, min(as.integer(ncores_grid), nrow(jobs)))
    if (verbose) {
      message(sprintf("Running %d (hyper, fold) jobs on %d cores.", nrow(jobs), ncores_grid))
    }

    job_results <- parallel::mclapply(
      seq_len(nrow(jobs)),
      function(ii) run_one_job(jobs[ii, , drop = FALSE]),
      mc.cores = ncores_grid
    )
  } else {
    job_results <- lapply(
      seq_len(nrow(jobs)),
      function(ii) run_one_job(jobs[ii, , drop = FALSE])
    )
  }

  cv_results <- do.call(rbind, job_results)

  ## --- summaries ---
  summary_fold <- stats::aggregate(
    cindex ~ fold + k + lambda + nu + lambdaW + lambdaH + alpha,
    data = cv_results,
    FUN  = function(x) mean(x, na.rm = TRUE)
  )
  names(summary_fold)[names(summary_fold) == "cindex"] <- "mean_cindex_fold"

  ## helper: standard error across folds
  se_fun <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) <= 1L) return(NA_real_)
    stats::sd(x) / sqrt(length(x))
  }

  ## Step 2a: mean across folds (of the fold-level init-averaged C-index)
  summary_mean <- stats::aggregate(
    mean_cindex_fold ~ k + lambda + nu + lambdaW + lambdaH + alpha,
    data = summary_fold,
    FUN  = function(x) mean(x, na.rm = TRUE)
  )
  names(summary_mean)[names(summary_mean) == "mean_cindex_fold"] <- "mean_cindex"

  ## Step 2b: SE across folds (of the fold-level init-averaged C-index)
  summary_se <- stats::aggregate(
    mean_cindex_fold ~ k + lambda + nu + lambdaW + lambdaH + alpha,
    data = summary_fold,
    FUN  = se_fun
  )
  names(summary_se)[names(summary_se) == "mean_cindex_fold"] <- "se_cindex"

  ## Combine mean and SE
  summary <- merge(
    summary_mean,
    summary_se,
    by   = c("k", "lambda", "nu", "lambdaW", "lambdaH", "alpha"),
    sort = FALSE
  )

  list(
    cv_results   = cv_results,   # raw: per fold, per init, per alpha, per hyper
    summary_fold = summary_fold, # step 1: mean across inits, per fold+alpha+hyper
    summary      = summary,      # step 2: mean across folds + SE, per alpha+hyper
    alpha_grid   = alpha_grid,
    hyper_grid   = hyper_grid,
    folds        = folds
  )
}
