#' @keywords internal
#' @title CV over hyperparameter grid with alpha warmstarts
#' @description
#' Internal helper for cross-validation over a grid of hyperparameters
#' \code{k, alpha, lambda, nu, lambdaW, lambdaH}, using warmstarts along
#' \code{alpha_grid} via \code{desurv_alpha_warmstart()} for each combination
#' of \code{k, lambda, nu, lambdaW, lambdaH}. Parallelism (if enabled) is
#' applied at the (hyperparameter, fold, init) job level; within each job,
#' \code{desurv_alpha_warmstart} is invoked for a single initialization so
#' nested parallelism is still avoided. Folds whose training subset has zero
#' observed events are skipped (emitting a warning) and contribute
#' \code{NA} C-indices.
#' @param dataset Optional vector identifying datasets/cohorts per sample,
#'   used to stratify folds jointly by event indicator and dataset.
#' @export
.desurv_cv_grid_warmstart <- function(
    X, y = NULL, d = NULL, dataset = NULL,
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
    ntop          = NULL,
    preprocess    = FALSE,
    ngene         = 1000,
    genes         = NULL,
    method_trans_train = c("rank", "quant", "none"),
    parallel_grid = FALSE,
    ncores_grid   = NULL,
    verbose       = TRUE
) {
  preprocess_flag <- isTRUE(preprocess)
  method_trans_train <- match.arg(method_trans_train)

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

  if (inherits(X, "desurv_data")) {
    if (preprocess_flag) {
      stop("preprocess = TRUE is not supported when `X` is already a 'desurv_data' object.")
    }
    X_full <- X$X
    y_full <- X$y
    d_full <- X$d
    dataset_full <- X$dataset
    if (is.null(dataset_full)) {
      dataset_full <- rep("_all", ncol(X_full))
    }
    if (!is.null(dataset)) {
      dataset_full <- as.character(dataset)
      if (length(dataset_full) != ncol(X_full)) {
        stop("`dataset` must have length equal to ncol(X).")
      }
    }
  } else {
    X_full <- as.matrix(X)
    y_full <- y
    d_full <- d
    tmp    <- .validate_desurv_data(
      X_full, y_full, d_full, k = hyper_grid$k[1L], dataset = dataset
    )
    X_full <- tmp$X; y_full <- tmp$y; d_full <- tmp$d
    dataset_full <- tmp$dataset
  }
  ## --- alpha grid ---
  alpha_grid <- as.numeric(alpha_grid)
  if (length(alpha_grid) == 0L || anyNA(alpha_grid)) {
    stop("alpha_grid must be a non-empty numeric vector without NA.")
  }
  A <- length(alpha_grid)

  ## --- n_starts / nfolds checks ---
  n_starts <- as.integer(n_starts)
  if (length(n_starts) != 1L || is.na(n_starts) || n_starts <= 0L) {
    stop("n_starts must be a positive integer.")
  }

  nfolds <- as.integer(nfolds)
  if (length(nfolds) != 1L || is.na(nfolds) || nfolds <= 1L) {
    stop("nfolds must be an integer >= 2.")
  }

  n <- ncol(X_full)

  ## --- folds: stratified by event ---
  if (is.null(folds)) {
    folds <- .desurv_make_folds_stratified(
      d_full,
      dataset = dataset_full,
      nfolds  = nfolds,
      seed    = seed
    )
    nfolds <- attr(folds, "nfolds")
    attr(folds, "nfolds") <- NULL
  } else {
    folds <- as.integer(folds)
    if (length(folds) != n) {
      stop("folds must have length equal to ncol(X).")
    }
    if (!all(folds %in% seq_len(nfolds))) {
      stop("All entries in folds must be in 1:nfolds.")
    }
    if (length(unique(folds)) != nfolds) {
      stop("Provided folds must contain all integers in 1:nfolds.")
    }
  }

  fold_levels <- sort(unique(folds))
  nfolds_effective <- length(fold_levels)

  if (verbose) {
    message(sprintf(
      "Running CV over %d hyperparameter combos, %d alphas, %d starts, %d folds.",
      H, A, n_starts, nfolds_effective
    ))
  }

  ## --- define (hyper, fold) jobs and cache payloads ---
  base_jobs <- expand.grid(
    h = seq_len(H),
    f = fold_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (nrow(base_jobs) == 0L) {
    stop("Hyperparameter grid is empty; specify at least one combination.")
  }
  base_jobs$job_id <- seq_len(nrow(base_jobs))

  fold_cache <- vector("list", length = length(fold_levels))
  names(fold_cache) <- as.character(fold_levels)
  for (f in fold_levels) {
    idx_val <- which(folds == f)
    idx_tr  <- setdiff(seq_len(n), idx_val)
    fold_cache[[as.character(f)]] <- .desurv_prepare_fold_payload(
      X_full             = X_full,
      y_full             = y_full,
      d_full             = d_full,
      dataset_full       = dataset_full,
      idx_tr             = idx_tr,
      idx_val            = idx_val,
      preprocess         = preprocess_flag,
      ngene              = ngene,
      genes              = genes,
      method_trans_train = method_trans_train,
      verbose            = verbose
    )
  }

  job_payloads <- vector("list", length = nrow(base_jobs))
  for (job_idx in seq_len(nrow(base_jobs))) {
    job_row <- base_jobs[job_idx, , drop = FALSE]
    h <- job_row$h
    f <- job_row$f

    hp <- hyper_grid[h, ]
    k_cur       <- hp$k
    lambda_cur  <- hp$lambda
    nu_cur      <- hp$nu
    lambdaW_cur <- hp$lambdaW
    lambdaH_cur <- hp$lambdaH

    fold_payload <- fold_cache[[as.character(f)]]
    has_events_tr <- sum(fold_payload$d_tr) > 0
    data_tr <- if (has_events_tr) {
      desurv_data(fold_payload$X_tr, fold_payload$y_tr, fold_payload$d_tr, k_cur)
    } else {
      NULL
    }

    seed_fold <- if (is.null(seed)) NULL else {
      as.integer(seed + 10000L * (h - 1L) + 100L * (f - 1L))
    }

    job_payloads[[job_idx]] <- list(
      k_cur       = k_cur,
      lambda_cur  = lambda_cur,
      nu_cur      = nu_cur,
      lambdaW_cur = lambdaW_cur,
      lambdaH_cur = lambdaH_cur,
      data_tr     = data_tr,
      X_val       = fold_payload$X_val,
      y_val       = fold_payload$y_val,
      d_val       = fold_payload$d_val,
      seed_fold   = seed_fold,
      has_events_tr = has_events_tr
    )
  }

  jobs <- base_jobs[rep(seq_len(nrow(base_jobs)), each = n_starts), , drop = FALSE]
  jobs$init_id <- rep(seq_len(n_starts), times = nrow(base_jobs))

  ## --- inner job function ---
  run_one_job <- function(job_row) {
    h <- job_row$h
    f <- job_row$f
    init_id <- job_row$init_id
    payload <- job_payloads[[job_row$job_id]]

    k_cur       <- payload$k_cur
    lambda_cur  <- payload$lambda_cur
    nu_cur      <- payload$nu_cur
    lambdaW_cur <- payload$lambdaW_cur
    lambdaH_cur <- payload$lambdaH_cur

    if (verbose && !parallel_grid) {
      if (init_id == 1L) {
        fold_pos <- match(f, fold_levels)
        message(sprintf(
          "Hyper combo %d/%d (k=%d, lambda=%.3g, nu=%.3g, lambdaW=%.3g, lambdaH=%.3g), fold %d/%d",
          h, H, k_cur, lambda_cur, nu_cur, lambdaW_cur, lambdaH_cur, fold_pos, nfolds_effective
        ))
      }
      message(sprintf("  Initialization %d/%d", init_id, n_starts))
    }

    if (!payload$has_events_tr) {
      if (init_id == 1L) {
        warning(sprintf(
          "Skipping hyper combo %d/%d, fold %d: training data has zero events; returning NA C-index.",
          h, H, f
        ), call. = FALSE)
      }
      df <- data.frame(
        fold    = f,
        k       = k_cur,
        lambda  = lambda_cur,
        nu      = nu_cur,
        lambdaW = lambdaW_cur,
        lambdaH = lambdaH_cur,
        init_id = init_id,
        alpha   = alpha_grid,
        cindex  = NA_real_,
        row.names = NULL
      )
      df <- df[, c("fold", "k", "lambda", "nu", "lambdaW", "lambdaH",
                   "init_id", "alpha", "cindex")]
      diag_df <- df
      diag_df$train_cindex <- NA_real_
      diag_df$train_convergence <- NA
      diag_df$train_nan_flag <- NA
      diag_df$val_has_events <- sum(payload$d_val) > 0L
      diag_df$val_lp_has_na <- NA
      diag_df$val_lp_has_inf <- NA
      diag_df$val_cindex <- NA_real_
      diag_df$status <- "no_train_events"
      diag_df$message <- "Training data has zero events"
      return(list(cv = df, diag = diag_df))
    }

    # unique seed for this (hyper, fold)
    seed_fold <- payload$seed_fold

    # warmstart path on TRAINING data for this combo (sequential inside)
    ws_res <- desurv_alpha_warmstart(
      X          = payload$data_tr,
      alpha_grid = alpha_grid,
      lambda     = lambda_cur,
      nu         = nu_cur,
      lambdaW    = lambdaW_cur,
      lambdaH    = lambdaH_cur,
      n_starts   = 1L,
      seed       = if (is.null(seed_fold)) NULL else as.integer(seed_fold + (init_id - 1L)),
      tol        = tol,
      maxit      = maxit,
      parallel   = FALSE,    # important: avoid nested parallelism
      ncores     = NULL,
      verbose    = verbose && !parallel_grid
    )

    # compute validation C-index vector for this init: length A
    cindex_vals <- rep(NA_real_, A)
    diag_rows <- vector("list", A)
    fits_one <- ws_res$fits[[1L]]
    val_has_events <- sum(payload$d_val) > 0L
    surv_obj_val <- survival::Surv(payload$y_val, payload$d_val)
    for (i in seq_len(A)) {
      fit_ji <- fits_one[[i]]
      lp_val <- .desurv_lp_with_top_genes(
        fit_ji,
        X_new = payload$X_val,
        ntop = ntop
      )

      lp_has_na <- anyNA(lp_val)
      lp_has_inf <- any(is.infinite(lp_val))
      lp_valid <- !(lp_has_na || lp_has_inf)

      diag_row <- data.frame(
        fold    = f,
        k       = k_cur,
        lambda  = lambda_cur,
        nu      = nu_cur,
        lambdaW = lambdaW_cur,
        lambdaH = lambdaH_cur,
        init_id = init_id,
        alpha   = alpha_grid[i],
        train_cindex = fit_ji$cindex,
        train_convergence = isTRUE(fit_ji$convergence),
        train_nan_flag = isTRUE(fit_ji$nan_flag),
        val_has_events = val_has_events,
        val_lp_has_na = lp_has_na,
        val_lp_has_inf = lp_has_inf,
        val_cindex = NA_real_,
        status = "ok",
        message = NA_character_,
        row.names = NULL,
        stringsAsFactors = FALSE
      )

      if (!val_has_events) {
        cindex_vals[i] <- NA_real_
        diag_row$status <- "no_val_events"
        diag_row$message <- "Validation fold has zero events"
      } else if (!lp_valid) {
        cindex_vals[i] <- NA_real_
        diag_row$status <- "nonfinite_linear_predictor"
        diag_row$message <- "Linear predictor contained non-finite values"
      } else {
        cindex_val <- cvwrapr::getCindex(lp_val, surv_obj_val)
        cindex_vals[i] <- cindex_val
        diag_row$val_cindex <- cindex_val
        if (is.na(cindex_val)) {
          diag_row$status <- "cindex_na"
          diag_row$message <- "C-index computation returned NA"
        }
      }
      diag_rows[[i]] <- diag_row
    }
    diag_df <- do.call(rbind, diag_rows)

    # tidy data.frame for this (hyper, fold, init) job
    df <- data.frame(
      fold    = f,
      k       = k_cur,
      lambda  = lambda_cur,
      nu      = nu_cur,
      lambdaW = lambdaW_cur,
      lambdaH = lambdaH_cur,
      init_id = init_id,
      alpha   = alpha_grid,
      cindex  = cindex_vals,
      row.names = NULL
    )

    # reorder columns a bit
    df <- df[, c("fold", "k", "lambda", "nu", "lambdaW", "lambdaH",
                 "init_id", "alpha", "cindex")]
    list(cv = df, diag = diag_df)
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
      message(sprintf(
        "Running %d (hyper, fold, init) jobs on %d cores.",
        nrow(jobs), ncores_grid
      ))
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

  cv_results <- do.call(rbind, lapply(job_results, `[[`, "cv"))
  cv_diagnostics <- do.call(rbind, lapply(job_results, `[[`, "diag"))

  ## --- summaries ---
  mean_fold_safe <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    mean(x)
  }

  summary_fold <- stats::aggregate(
    cindex ~ fold + k + lambda + nu + lambdaW + lambdaH + alpha,
    data = cv_results,
    FUN  = mean_fold_safe,
    na.action = stats::na.pass
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
    FUN  = function(x) {
      if (anyNA(x) || any(!is.finite(x))) return(NA_real_)
      mean(x)
    },
    na.action = stats::na.pass
  )
  names(summary_mean)[names(summary_mean) == "mean_cindex_fold"] <- "mean_cindex"

  ## Step 2b: SE across folds (of the fold-level init-averaged C-index)
  summary_se <- stats::aggregate(
    mean_cindex_fold ~ k + lambda + nu + lambdaW + lambdaH + alpha,
    data = summary_fold,
    FUN  = se_fun,
    na.action = stats::na.pass
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
    folds        = folds,
    diagnostics  = cv_diagnostics
  )
}
