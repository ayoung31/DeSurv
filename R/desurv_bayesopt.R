#' Bayesian optimisation of DeSurv hyperparameters
#'
#' @description
#' Sequentially tunes the scalar inputs to [desurv_cv()] with a simple Gaussian
#' process surrogate and expected-improvement acquisition rule. Each evaluation
#' runs `desurv_cv(cv_only = TRUE)` for a single hyperparameter combination, so
#' the optimisation surface is defined by the resulting mean validation
#' C-index. The routine always maximises the C-index.
#'
#' @param X,y,d,dataset,samp_keeps Inputs forwarded to [desurv_cv()]. Provide
#'   `dataset`/`samp_keeps` only when you plan to enable preprocessing inside
#'   cross-validation.
#' @param preprocess,method_trans_train,engine Arguments passed through to
#'   [desurv_cv()]. Bayesian optimisation assumes the output metric is
#'   comparable across iterations, so keep these unchanged throughout the run.
#' @param bo_bounds Named list describing the search space. Each entry must be
#'   either a numeric vector of length two (`c(lower, upper)`) or a list with
#'   elements `lower`, `upper`, and optional `type = "continuous"`/`"integer"`
#'   plus `scale = "linear"`/`"log10"`. Defaults draw from the ranges supplied
#'   in the prompt (see [desurv_bo_default_bounds()]).
#' @param bo_fixed Named list of scalar arguments that should remain constant
#'   during optimisation (e.g., `lambdaW_grid = 0`, `lambdaH_grid = 0`,
#'   `nfolds = 5`). These override the defaults used in the optimiser.
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
#' @param ntop Optional positive integer passed to [desurv_cv()] to restrict
#'   each factor to its top genes when scoring validation folds. Set to
#'   `NULL` or a non-positive value to disable the filtering.
#' @param verbose Logical; print progress after each evaluation/iteration.
#' @param ... Additional arguments forwarded to [desurv_cv()] (for example
#'   `seed`, `folds`, or `parallel_grid`).
#'
#' @details
#' This helper relies on the `DiceKriging` and `lhs` packages for Gaussian
#' process modelling and Latin hypercube sampling. They are listed under
#' `Suggests`; an informative error is thrown if they are not installed.
#'
#' @return An object of class `desurv_bo` containing:
#' \itemize{
#'   \item `history`: data frame of every evaluation (parameters, score, status).
#'   \item `best`: list with `params` (named numeric vector) and
#'         `mean_cindex`.
#'   \item `bounds`: parsed search-space specification.
#'   \item `fixed`: fixed hyperparameters supplied through `bo_fixed`.
#'   \item `seed`: RNG seed used internally.
#'   \item `call`: original function call.
#'   \item `km_fit`: last successfully fitted Gaussian-process surrogate (or `NULL` if skipped).
#' }
#'
#' @seealso [desurv_cv()], [desurv_bo_default_bounds()]
#'
#' @examples
#' \dontrun{
#' data <- readRDS("data.rds")
#' bo_fit <- desurv_bayesopt(
#'   X = data$ex,
#'   y = data$sampInfo$time,
#'   d = data$sampInfo$event,
#'   dataset = data$sampInfo$dataset,
#'   samp_keeps = data$samp_keeps,
#'   preprocess = TRUE,
#'   method_trans_train = "rank",
#'   engine = "cold",
#'   bo_bounds = desurv_bo_default_bounds(include_factor_penalties = FALSE),
#'   bo_fixed = list(lambdaW_grid = 0, lambdaH_grid = 0, nfolds = 5),
#'   n_init = 12,
#'   n_iter = 25,
#'   candidate_pool = 2000,
#'   exploration_weight = 0.01,
#'   seed = 2025,
#'   verbose = TRUE
#' )
#' print(bo_fit)
#' }
#'
#' @export
desurv_bayesopt <- function(
    X,
    y,
    d,
    dataset = NULL,
    samp_keeps = NULL,
    preprocess = TRUE,
    method_trans_train = c("rank", "quant", "none"),
    engine = c("cold", "warmstart"),
    bo_bounds = desurv_bo_default_bounds(),
    bo_fixed = list(),
    n_init = NULL,
    n_iter = 20,
    candidate_pool = NULL,
    exploration_weight = 0,
    seed = NULL,
    cv_verbose = FALSE,
    ntop = NULL,
    verbose = TRUE,
    ...
) {
  method_trans_train <- match.arg(method_trans_train)
  engine <- match.arg(engine)

  if (!requireNamespace("DiceKriging", quietly = TRUE)) {
    stop("Package 'DiceKriging' is required for desurv_bayesopt(). Please install it.",
         call. = FALSE)
  }
  if (!requireNamespace("lhs", quietly = TRUE)) {
    stop("Package 'lhs' is required for desurv_bayesopt(). Please install it.",
         call. = FALSE)
  }

  bound_info <- .desurv_bo_parse_bounds(bo_bounds)
  if (nrow(bound_info) == 0L) {
    stop("`bo_bounds` must contain at least one tunable parameter.", call. = FALSE)
  }

  d_dim <- nrow(bound_info)
  n_init <- n_init %||% max(5L, 3L * d_dim)
  candidate_pool <- candidate_pool %||% max(1000L, 200L * d_dim)
  exploration_weight <- max(0, as.numeric(exploration_weight))

  rng_seed <- if (is.null(seed)) sample.int(1e6, 1L) else as.integer(seed)
  set.seed(rng_seed)

  base_args <- list(
    X = X,
    y = y,
    d = d,
    dataset = dataset,
    samp_keeps = samp_keeps,
    preprocess = preprocess,
    method_trans_train = method_trans_train,
    engine = engine,
    cv_only = TRUE,
    ntop = ntop,
    verbose = isTRUE(cv_verbose),
    ...
  )

  defaults_fixed <- list(
    lambdaW_grid = 0,
    lambdaH_grid = 0,
    nfolds = 5,
    n_starts = 20,
    ngene = if (preprocess) 2000 else NULL,
    tol = 1e-4,
    maxit = 1000
  )
  defaults_fixed <- defaults_fixed[!vapply(defaults_fixed, is.null, logical(1))]
  fixed_args <- utils::modifyList(defaults_fixed, bo_fixed, keep.null = TRUE)

  existing_keys <- character()
  history_rows <- list()
  unit_store <- list()
  eval_id <- 0L
  last_km_fit <- NULL

  eval_point <- function(point, stage, iter) {
    eval_id <<- eval_id + 1L
    start_time <- proc.time()[3]

    args <- .desurv_merge_args(base_args, fixed_args, point$values)
    result <- tryCatch(do.call(desurv_cv, args), error = function(e) e)
    elapsed <- proc.time()[3] - start_time

    if (inherits(result, "error")) {
      score <- NA_real_
      status <- "error"
      msg <- conditionMessage(result)
    } else {
      metric <- result$summary$mean_cindex
      score <- metric[1L]
      status <- "ok"
      msg <- NA_character_
    }

    unit_store[[eval_id]] <<- point$unit

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

  lhs_init <- lhs::randomLHS(n_init, d_dim)
  for (i in seq_len(n_init)) {
    pt <- draw_point(lhs_init[i, ])
    if (is.null(pt)) next
    eval_point(pt, stage = "initial", iter = 0L)
  }

  if (!length(history_rows)) {
    stop("No evaluations were completed. Check `bo_bounds` and the supplied data.", call. = FALSE)
  }

  for (iter in seq_len(n_iter)) {
    hist_df <- do.call(rbind, history_rows)
    ok_idx <- hist_df$status == "ok" & !is.na(hist_df$mean_cindex)
    if (sum(ok_idx) < 2L) {
      pt <- draw_point()
      if (is.null(pt)) break
      eval_point(pt, stage = "random", iter = iter)
      next
    }

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
      error = function(e) NULL
    )
    last_km_fit <- km_fit

    candidate_units <- lhs::randomLHS(candidate_pool, d_dim)
    colnames(candidate_units) <- bound_info$parameter

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
        error = function(e) NULL
      )
      if (is.null(preds)) {
        return(draw_point())
      }

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

      order_idx <- order(ei, decreasing = TRUE)
      for (idx in order_idx) {
        pt <- draw_point(candidate_units[idx, ])
        if (!is.null(pt)) return(pt)
      }
      draw_point()
    }

    next_pt <- acquire_point()
    if (is.null(next_pt)) break
    eval_point(next_pt, stage = "bo", iter = iter)
  }

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

  structure(
    list(
      history = history_df,
      best = list(
        params = best_params,
        mean_cindex = best_row$mean_cindex
      ),
      bounds = bound_info,
      fixed = c(fixed_args, list(ntop = ntop)),
      seed = rng_seed,
      call = match.call(),
      km_fit = last_km_fit
    ),
    class = "desurv_bo"
  )
}

#' Default Bayesian optimisation bounds for DeSurv
#'
#' @description
#' Convenience helper returning the recommended tuning ranges presented in the
#' user prompt. Pass `include_factor_penalties = FALSE` to keep `lambdaW_grid`
#' and `lambdaH_grid` fixed (e.g., at zero via `bo_fixed`).
#'
#' @param include_factor_penalties Logical; include search ranges for
#'   `lambdaW_grid`/`lambdaH_grid`. Defaults to `TRUE`.
#'
#' @return Named list suitable for the `bo_bounds` argument of
#'   [desurv_bayesopt()].
#'
#' @export
desurv_bo_default_bounds <- function(include_factor_penalties = TRUE) {
  bounds <- list(
    k_grid = list(lower = 2, upper = 12, type = "integer"),
    alpha_grid = list(lower = 0, upper = 0.95, type = "continuous"),
    lambda_grid = list(lower = 1e-5, upper = 1e5, scale = "log10"),
    nu_grid = list(lower = 0, upper = 1),
    n_starts = list(lower = 10, upper = 100, type = "integer"),
    nfolds = list(lower = 5, upper = 10, type = "integer"),
    ngene = list(lower = 1000, upper = 5000, type = "integer"),
    tol = list(lower = 1e-6, upper = 1e-4, scale = "log10"),
    maxit = list(lower = 200, upper = 6000, type = "integer"),
    ntop = list(lower = 25, upper = 100, type = "integer")
  )
  if (isTRUE(include_factor_penalties)) {
    bounds$lambdaW_grid <- list(lower = 1e-5, upper = 1e5, scale = "log10")
    bounds$lambdaH_grid <- list(lower = 1e-5, upper = 1e5, scale = "log10")
  }
  bounds
}

#' @export
print.desurv_bo <- function(x, ...) {
  cat("Bayesian optimisation for DeSurv hyperparameters\n")
  cat("Evaluations :", nrow(x$history), "\n")
  cat("Best C-index:", sprintf("%.4f", x$best$mean_cindex), "\n")
  cat("Best params :\n")
  print(x$best$params)
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.desurv_bo_parse_bounds <- function(bounds) {
  if (!length(bounds)) {
    return(data.frame())
  }
  param <- names(bounds)
  info <- lapply(seq_along(bounds), function(i) {
    entry <- bounds[[i]]
    if (is.numeric(entry) && length(entry) == 2L) {
      list(lower = entry[1], upper = entry[2], type = "continuous", scale = "linear")
    } else if (is.list(entry)) {
      lower <- entry$lower
      upper <- entry$upper
      type <- entry$type %||% "continuous"
      scale <- entry$scale %||% "linear"
      list(lower = lower, upper = upper, type = type, scale = scale)
    } else {
      stop("Each element of `bo_bounds` must be a numeric length-2 vector or a list.",
           call. = FALSE)
    }
  })
  lower <- vapply(info, function(x) x$lower, numeric(1))
  upper <- vapply(info, function(x) x$upper, numeric(1))
  if (any(!is.finite(lower) | !is.finite(upper))) {
    stop("Lower/upper bounds must be finite.", call. = FALSE)
  }
  if (any(upper < lower)) {
    stop("Each lower bound must be <= the corresponding upper bound.", call. = FALSE)
  }
  type <- vapply(info, function(x) match.arg(x$type, c("continuous", "integer")), character(1))
  scale <- vapply(info, function(x) match.arg(x$scale, c("linear", "log10")), character(1))
  if (any(scale == "log10" & lower <= 0)) {
    stop("log10-scaled bounds must be strictly positive.", call. = FALSE)
  }
  data.frame(
    parameter = param,
    lower = lower,
    upper = upper,
    type = type,
    scale = scale,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.desurv_bo_make_point <- function(unit_vec, info) {
  values <- mapply(
    function(u, lower, upper, type, scale) {
      if (lower == upper) {
        val <- lower
      } else if (scale == "log10") {
        log_lower <- log10(lower)
        log_upper <- log10(upper)
        val <- 10^(log_lower + u * (log_upper - log_lower))
      } else {
        val <- lower + u * (upper - lower)
      }
      if (type == "integer") {
        val <- round(val)
      }
      val <- max(lower, min(upper, val))
      val
    },
    u = as.list(unit_vec),
    lower = info$lower,
    upper = info$upper,
    type = info$type,
    scale = info$scale,
    SIMPLIFY = TRUE
  )
  names(values) <- info$parameter
  list(
    unit = unit_vec,
    values = as.list(values),
    key = paste(sprintf("%s=%s", names(values), signif(values, 8)), collapse = "|")
  )
}
