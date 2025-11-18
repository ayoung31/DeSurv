#' Bayesian optimisation of DeSurv CV hyperparameters
#'
#' @description
#' Sequentially tunes hyperparameters for [desurv_cv()] using Bayesian
#' optimization with a Gaussian process surrogate and expected-improvement
#' acquisition rule. Unlike grid search which evaluates all combinations,
#' this function efficiently explores the hyperparameter space by calling
#' `desurv_cv()` with single parameter values (not grids) and maximizing
#' the resulting mean validation C-index.
#'
#' @param X,y,d,dataset,samp_keeps Inputs forwarded to [desurv_cv()]. Provide
#'   `dataset`/`samp_keeps` when preprocessing is enabled.
#' @param preprocess,method_trans_train,engine Arguments passed through to
#'   [desurv_cv()]. Bayesian optimisation assumes the output metric is
#'   comparable across iterations, so keep these unchanged throughout the run.
#' @param nfolds,tol,maxit Fixed CV/optimization parameters passed to
#'   [desurv_cv()]. These can also be tuned via `bo_bounds` if desired, but
#'   typically you want to fix them for computational efficiency.
#' @param ntop Optional positive integer passed to [desurv_cv()] to restrict
#'   the validation predictor to the top genes per factor. Set to `NULL` or a
#'   non-positive value to disable this filtering.
#' @param bo_bounds Named list describing the search space for hyperparameters
#'   to tune. Each entry must be either a numeric vector of length two
#'   (`c(lower, upper)`) or a list with elements `lower`, `upper`, and
#'   optional `type = "continuous"`/`"integer"` plus
#'   `scale = "linear"`/`"log10"`. See [desurv_cv_bo_default_bounds()] for
#'   recommended ranges.
#' @param n_init Number of Latin-hypercube evaluations before fitting the GP.
#'   Defaults to `max(5, 3 * length(bo_bounds))`.
#' @param n_iter Number of sequential BO iterations after the initial design.
#' @param candidate_pool Number of random candidates assessed per acquisition
#'   step. Defaults to `max(1000, 200 * length(bo_bounds))`.
#' @param exploration_weight Non-negative scalar `xi` used in the expected
#'   improvement formula to encourage exploration.
#' @param seed Optional base seed for reproducibility.
#' @param cv_verbose Logical; if `TRUE`, propagate `verbose = TRUE` into
#'   `desurv_cv()` calls. The optimiser itself still reports progress via
#'   `verbose`.
#' @param verbose Logical; print progress after each evaluation/iteration.
#' @param ... Additional arguments forwarded to [desurv_cv()] (for example
#'   `folds`, `parallel_grid`, `ncores_grid`).
#'
#' @details
#' This function optimizes hyperparameters like `k`, `alpha`, `lambda`, `nu`,
#' `lambdaW`, `lambdaH`, `n_starts`, and `ngene` by calling `desurv_cv()`
#' with single values (not grids) for each parameter. Each BO evaluation runs
#' cross-validation and extracts the mean C-index across folds.
#'
#' The helper relies on the `DiceKriging` and `lhs` packages for Gaussian
#' process modelling and Latin hypercube sampling. They are listed under
#' `Suggests`; an informative error is thrown if they are not installed.
#'
#' **Important**: The parameter names in `bo_bounds` should match the
#' arguments to `desurv_cv()` (e.g., `k_grid`, `alpha_grid`, etc.) even
#' though you're passing single values. The function will extract these
#' and pass them as scalars.
#'
#' @return An object of class `desurv_cv_bo` containing:
#' \itemize{
#'   \item `history`: data frame of every evaluation (parameters, score, status).
#'   \item `best`: list with `params` (named numeric vector) and
#'         `mean_cindex`.
#'   \item `bounds`: parsed search-space specification.
#'   \item `fixed`: fixed hyperparameters (nfolds, tol, maxit).
#'   \item `seed`: RNG seed used internally.
#'   \item `call`: original function call.
#'   \item `km_fit`: last successfully fitted Gaussian-process surrogate (or `NULL` if unavailable).
#' }
#'
#' @seealso [desurv_cv()], [desurv_cv_bo_default_bounds()]
#'
#' @examples
#' \dontrun{
#' data <- readRDS("data.rds")
#' bo_result <- desurv_cv_bayesopt(
#'   X = data$ex,
#'   y = data$sampInfo$time,
#'   d = data$sampInfo$event,
#'   dataset = data$sampInfo$dataset,
#'   samp_keeps = data$samp_keeps,
#'   preprocess = TRUE,
#'   method_trans_train = "rank",
#'   engine = "cold",
#'   nfolds = 5,
#'   tol = 1e-4,
#'   maxit = 100,
#'   bo_bounds = desurv_cv_bo_default_bounds(),
#'   n_init = 15,
#'   n_iter = 20,
#'   exploration_weight = 0.01,
#'   seed = 2025,
#'   verbose = TRUE
#' )
#' print(bo_result)
#'
#' # Extract best parameters
#' best_params <- bo_result$best$params
#'
#' # Refit with best parameters on full data
#' final_cv <- desurv_cv(
#'   X = data$ex,
#'   y = data$sampInfo$time,
#'   d = data$sampInfo$event,
#'   dataset = data$sampInfo$dataset,
#'   samp_keeps = data$samp_keeps,
#'   preprocess = TRUE,
#'   method_trans_train = "rank",
#'   k_grid = best_params["k_grid"],
#'   alpha_grid = best_params["alpha_grid"],
#'   lambda_grid = best_params["lambda_grid"],
#'   nu_grid = best_params["nu_grid"],
#'   lambdaW_grid = best_params["lambdaW_grid"],
#'   lambdaH_grid = best_params["lambdaH_grid"],
#'   n_starts = best_params["n_starts"],
#'   ngene = best_params["ngene"],
#'   nfolds = 5,
#'   tol = 1e-4,
#'   maxit = 100,
#'   engine = "cold"
#' )
#' }
#'
#' @export
desurv_cv_bayesopt <- function(
    X,
    y,
    d,
    dataset = NULL,
    samp_keeps = NULL,
    preprocess = TRUE,
    method_trans_train = c("rank", "quant", "none"),
    engine = c("cold", "warmstart"),
    nfolds = 5,
    tol = 1e-4,
    maxit = 100,
    ntop = NULL,
    bo_bounds = desurv_cv_bo_default_bounds(),
    n_init = NULL,
    n_iter = 20,
    candidate_pool = NULL,
    exploration_weight = 0,
    seed = NULL,
    cv_verbose = FALSE,
    verbose = TRUE,
    ...
) {
  method_trans_train <- match.arg(method_trans_train)
  engine <- match.arg(engine)

  # Check required packages
  if (!requireNamespace("DiceKriging", quietly = TRUE)) {
    stop("Package 'DiceKriging' is required for desurv_cv_bayesopt(). Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("Package 'lhs' is required for desurv_cv_bayesopt(). Please install it.",
         call. = FALSE)
  }

  # Parse and validate bounds
  bound_info <- .desurv_bo_parse_bounds(bo_bounds)
  if (nrow(bound_info) == 0L) {
    stop("`bo_bounds` must contain at least one tunable parameter.", call. = FALSE)
  }

  # Set defaults based on dimensionality
  d_dim <- nrow(bound_info)
  n_init <- n_init %||% max(5L, 3L * d_dim)
  candidate_pool <- candidate_pool %||% max(1000L, 200L * d_dim)
  exploration_weight <- max(0, as.numeric(exploration_weight))

  # Initialize RNG
  rng_seed <- if (is.null(seed)) sample.int(1e6, 1L) else as.integer(seed)
  set.seed(rng_seed)

  # Base arguments passed to every desurv_cv() call
  base_args <- list(
    X = X,
    y = y,
    d = d,
    dataset = dataset,
    samp_keeps = samp_keeps,
    preprocess = preprocess,
    method_trans_train = method_trans_train,
    engine = engine,
    nfolds = nfolds,
    tol = tol,
    maxit = maxit,
    ntop = ntop,
    cv_only = TRUE,  # Critical: only run CV, don't refit final model
    verbose = isTRUE(cv_verbose),
    ...
  )

  # Store fixed parameters for reporting
  fixed_params <- list(
    nfolds = nfolds,
    tol = tol,
    maxit = maxit,
    engine = engine,
    method_trans_train = method_trans_train,
    preprocess = preprocess,
    ntop = ntop
  )

  # Tracking variables
  existing_keys <- character()
  history_rows <- list()
  unit_store <- list()
  eval_id <- 0L
  last_km_fit <- NULL

  # Objective function: evaluate a single parameter configuration
  eval_point <- function(point, stage, iter) {
    eval_id <<- eval_id + 1L
    start_time <- proc.time()[3]

    # Combine base args with current parameter values
    args <- .desurv_merge_args(base_args, point$values)
    result <- tryCatch(do.call(desurv_cv, args), error = function(e) e)
    elapsed <- proc.time()[3] - start_time

    if (inherits(result, "error")) {
      score <- NA_real_
      status <- "error"
      msg <- conditionMessage(result)
    } else {
      # Extract mean C-index from CV results
      metric <- result$summary$mean_cindex
      score <- metric[1L]
      status <- "ok"
      msg <- NA_character_
    }

    # Store unit coordinates for GP fitting
    unit_store[[eval_id]] <<- point$unit

    # Record evaluation in history
    param_df <- as.data.frame(as.list(point$values), optional = TRUE)
    row <- cbind(
      data.frame(
        eval_id = eval_id,
        stage = stage,
        iteration = iter,
        stringsAsFactors = FALSE
      ),
      param_df,
      data.frame(
        mean_cindex = score,
        status = status,
        message = msg,
        elapsed = elapsed,
        stringsAsFactors = FALSE
      )
    )
    history_rows[[length(history_rows) + 1L]] <<- row

    # Print progress
    if (verbose) {
      msg_core <- sprintf(
        "[%s] eval %d (iter %d) -> mean_cindex = %s",
        stage,
        eval_id,
        iter,
        if (is.na(score)) "NA" else sprintf("%.4f", score)
      )
      if (!is.na(msg)) {
        message(msg_core, " | error: ", msg)
      } else {
        message(msg_core)
      }
    }

    list(score = score, status = status)
  }

  # Sample a point with collision detection
  draw_point <- function(unit_vec = NULL, max_attempts = 100L) {
    for (attempt in seq_len(max_attempts)) {
      u <- if (is.null(unit_vec)) stats::runif(d_dim) else unit_vec
      point <- .desurv_bo_make_point(u, bound_info)
      if (!point$key %in% existing_keys) {
        existing_keys <<- c(existing_keys, point$key)
        return(point)
      }
      unit_vec <- NULL
    }
    NULL
  }

  # Phase 1: Latin Hypercube Sampling initialization
  if (verbose) message("Starting BO initialization with ", n_init, " LHS samples...")
  lhs_init <- lhs::randomLHS(n_init, d_dim)
  for (i in seq_len(n_init)) {
    pt <- draw_point(lhs_init[i, ])
    if (is.null(pt)) next
    eval_point(pt, stage = "initial", iter = 0L)
  }

  if (!length(history_rows)) {
    stop("No evaluations were completed. Check `bo_bounds` and the supplied data.", call. = FALSE)
  }

  # Phase 2: Sequential BO iterations
  if (verbose) message("Starting ", n_iter, " BO iterations...")
  for (iter in seq_len(n_iter)) {
    hist_df <- do.call(rbind, history_rows)
    ok_idx <- hist_df$status == "ok" & !is.na(hist_df$mean_cindex)

    # Need at least 2 points to fit GP
    if (sum(ok_idx) < 2L) {
      if (verbose) message("  Iter ", iter, ": insufficient data for GP, sampling randomly")
      pt <- draw_point()
      if (is.null(pt)) break
      eval_point(pt, stage = "random", iter = iter)
      next
    }

    # Fit Gaussian Process surrogate
    design <- do.call(rbind, unit_store[hist_df$eval_id[ok_idx]])
    colnames(design) <- bound_info$parameter
    response <- hist_df$mean_cindex[ok_idx]

    km_fit <- tryCatch(
      DiceKriging::km(
        design = design,
        response = response,
        covtype = "matern5_2",
        nugget.estim = TRUE,
        control = list(trace = FALSE)
      ),
      error = function(e) {
        if (verbose) message("  Iter ", iter, ": GP fit failed, will sample randomly")
        NULL
      }
    )
    last_km_fit <- km_fit

    # Generate candidate pool for acquisition function
    candidate_units <- lhs::randomLHS(candidate_pool, d_dim)
    colnames(candidate_units) <- bound_info$parameter

    # Acquisition function: Expected Improvement
    acquire_point <- function() {
      if (is.null(km_fit)) {
        return(draw_point())
      }

      preds <- tryCatch(
        DiceKriging::predict.km(
          km_fit,
          newdata = candidate_units,
          type = "UK",
          checkNames = FALSE,
          se.compute = TRUE,
          cov.compute = FALSE
        ),
        error = function(e) {
          if (verbose) message("  Iter ", iter, ": GP prediction failed")
          NULL
        }
      )
      if (is.null(preds)) {
        return(draw_point())
      }

      # Expected Improvement calculation
      current_best <- max(response)
      improv <- preds$mean - current_best - exploration_weight
      sd <- preds$sd
      with_sd <- sd > 0
      z <- rep(0, length(sd))
      z[with_sd] <- improv[with_sd] / sd[with_sd]
      ei <- numeric(length(sd))
      ei[with_sd] <- improv[with_sd] * stats::pnorm(z[with_sd]) +
        sd[with_sd] * stats::dnorm(z[with_sd])
      ei[!with_sd] <- 0

      # Select candidate with highest EI (with collision avoidance)
      order_idx <- order(ei, decreasing = TRUE)
      for (idx in order_idx) {
        pt <- draw_point(candidate_units[idx, ])
        if (!is.null(pt)) return(pt)
      }
      draw_point()
    }

    next_pt <- acquire_point()
    if (is.null(next_pt)) {
      if (verbose) message("  Iter ", iter, ": no valid candidates, stopping early")
      break
    }
    eval_point(next_pt, stage = "bo", iter = iter)
  }

  # Extract best result
  history_df <- do.call(rbind, history_rows)
  ok <- history_df$status == "ok" & !is.na(history_df$mean_cindex)
  if (!any(ok)) {
    stop("Bayesian optimisation failed: no successful desurv_cv() evaluations.", call. = FALSE)
  }

  best_idx <- which.max(history_df$mean_cindex[ok])
  ok_rows <- history_df[ok, , drop = FALSE]
  best_row <- ok_rows[best_idx, , drop = FALSE]
  param_names <- bound_info$parameter
  best_params <- as.numeric(best_row[1, param_names, drop = TRUE])
  names(best_params) <- param_names

  if (verbose) {
    message("\nBayesian optimization complete!")
    message("Best mean C-index: ", sprintf("%.4f", best_row$mean_cindex))
    message("Best parameters:")
    for (pname in param_names) {
      message("  ", pname, ": ", best_params[pname])
    }
  }

  structure(
    list(
      history = history_df,
      best = list(
        params = best_params,
        mean_cindex = best_row$mean_cindex
      ),
      bounds = bound_info,
      fixed = fixed_params,
      seed = rng_seed,
      call = match.call(),
      km_fit = last_km_fit
    ),
    class = "desurv_cv_bo"
  )
}

#' Default Bayesian optimisation bounds for desurv_cv hyperparameters
#'
#' @description
#' Convenience helper returning the recommended tuning ranges for `desurv_cv()`
#' hyperparameters. These ranges are based on empirical experience and the
#' grid search results provided in the user prompt.
#'
#' @param include_factor_penalties Logical; include search ranges for
#'   `lambdaW_grid`/`lambdaH_grid`. Defaults to `TRUE`.
#' @param include_n_starts Logical; include `n_starts` in the search space.
#'   Defaults to `TRUE`. Set to `FALSE` and fix via function argument if you
#'   want to control computational cost.
#' @param include_ngene Logical; include `ngene` (number of genes for
#'   preprocessing) in the search space. Defaults to `TRUE`. Only relevant
#'   when `preprocess = TRUE`.
#'
#' @return Named list suitable for the `bo_bounds` argument of
#'   [desurv_cv_bayesopt()].
#'
#' @export
desurv_cv_bo_default_bounds <- function(
    include_factor_penalties = TRUE,
    include_n_starts = TRUE,
    include_ngene = TRUE
) {
  bounds <- list(
    k_grid = list(lower = 2, upper = 12, type = "integer"),
    alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
    lambda_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
    nu_grid = list(lower = 0, upper = 1, type = "continuous"),
    ntop = list(lower = 25, upper = 100, type = "integer")
  )

  if (isTRUE(include_factor_penalties)) {
    bounds$lambdaW_grid <- list(lower = 1e-5, upper = 1e5, scale = "log10")
    bounds$lambdaH_grid <- list(lower = 1e-5, upper = 1e5, scale = "log10")
  }

  if (isTRUE(include_n_starts)) {
    bounds$n_starts <- list(lower = 10, upper = 100, type = "integer")
  }

  if (isTRUE(include_ngene)) {
    bounds$ngene <- list(lower = 1000, upper = 5000, type = "integer")
  }

  bounds
}

#' @export
print.desurv_cv_bo <- function(x, ...) {
  cat("Bayesian optimisation for desurv_cv hyperparameters\n")
  cat("Evaluations :", nrow(x$history), "\n")
  cat("Best C-index:", sprintf("%.4f", x$best$mean_cindex), "\n")
  cat("Best params :\n")
  print(x$best$params)
  invisible(x)
}
