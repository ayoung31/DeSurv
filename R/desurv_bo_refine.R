 #' Refine Bayesian optimisation bounds based on prior history
#'
#' @description
#' Shrinks an existing set of BO bounds using the top-performing evaluations
#' from a previous run while respecting the original search space. Parameter
#' importance scores can be supplied to tighten high-impact parameters more
#' aggressively than low-impact ones.
#'
#' @param bounds Named list describing the current search space (same format
#'   as `bo_bounds` passed to [desurv_bayesopt()] or [desurv_cv_bayesopt()]).
#' @param history Data frame of past evaluations containing the parameter
#'   columns plus `mean_cindex` and `status`.
#' @param top_k Integer; number of best evaluations to consider when defining
#'   the allowable refined range. Defaults to 10.
#' @param importance Optional named numeric vector of parameter importances
#'   scaled to `[0, 1]`. Names must match the parameter keys in `bounds`. When
#'   `NULL`, all parameters are treated equally.
#' @param shrink_base Baseline shrinkage strength between 0 and 1. Values
#'   closer to 1 push the bounds quickly towards the top-k envelope.
#' @param importance_gain Additional shrinkage applied to each parameter
#'   according to its importance score (0-1). For example, with
#'   `importance_gain = 0.5`, a parameter with `importance = 1` will shrink
#'   `50%` faster than the base rate.
#'
#' @return A list of refined bounds ready to be passed back into the BO helper.
#'
#' @export
desurv_bo_refine_bounds <- function(
    bounds,
    history,
    top_k = 10,
    importance = NULL,
    shrink_base = 0.5,
    importance_gain = 0.3
) {
  if (!length(bounds)) {
    stop("`bounds` must contain at least one parameter.", call. = FALSE)
  }
  hist_df <- as.data.frame(history, stringsAsFactors = FALSE)
  if (!nrow(hist_df)) {
    warning("History is empty; returning the original bounds.")
    return(bounds)
  }

  info <- .desurv_bo_parse_bounds(bounds)
  params <- info$parameter
  missing_cols <- setdiff(params, names(hist_df))
  if (length(missing_cols)) {
    stop("History is missing parameter columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  ok_idx <- hist_df$status == "ok" & !is.na(hist_df$mean_cindex)
  if (!any(ok_idx)) {
    warning("History contains no successful evaluations; returning the original bounds.")
    return(bounds)
  }

  hist_ok <- hist_df[ok_idx, c("mean_cindex", params), drop = FALSE]
  hist_ok <- hist_ok[order(hist_ok$mean_cindex, decreasing = TRUE), , drop = FALSE]
  if (nrow(hist_ok) > top_k) {
    hist_ok <- hist_ok[seq_len(top_k), , drop = FALSE]
  }

  imp <- rep(1, length(params))
  names(imp) <- params
  if (!is.null(importance)) {
    keep <- intersect(names(importance), params)
    if (!length(keep)) {
      warning("Supplied `importance` vector had no matching parameter names; falling back to equal weights.")
    } else {
      imp[keep] <- importance[keep]
    }
  }
  imp <- pmax(0, pmin(1, imp))
  shrink_base <- min(max(shrink_base, 0), 1)
  importance_gain <- max(0, importance_gain)

  refined <- vector("list", length(params))
  names(refined) <- params

  for (i in seq_along(params)) {
    param <- params[i]
    entry <- bounds[[param]]
    info_row <- info[i, , drop = FALSE]
    current <- .desurv_bo_normalise_bound_entry(entry, info_row)
    top_vals <- hist_ok[[param]]
    env_low <- max(info_row$lower, min(top_vals, na.rm = TRUE))
    env_high <- min(info_row$upper, max(top_vals, na.rm = TRUE))

    if (!is.finite(env_low) || !is.finite(env_high) || env_low > env_high) {
      env_low <- env_high <- hist_ok[[param]][which.max(hist_ok$mean_cindex)]
    }

    shrink_weight <- min(1, shrink_base + importance_gain * imp[i])
    new_lower <- current$lower + shrink_weight * (env_low - current$lower)
    new_upper <- current$upper + shrink_weight * (env_high - current$upper)

    # Enforce the top-k envelope before applying type-specific tweaks
    new_lower <- max(min(new_lower, env_high), env_low)
    new_upper <- max(min(new_upper, env_high), env_low)

    if (identical(info_row$type, "integer")) {
      new_lower <- round(new_lower)
      new_upper <- round(new_upper)
    }
    new_lower <- max(min(new_lower, env_high), env_low)
    new_upper <- max(min(new_upper, env_high), env_low)
    new_lower <- max(info_row$lower, min(new_lower, info_row$upper))
    new_upper <- max(info_row$lower, min(new_upper, info_row$upper))
    if (new_lower > new_upper) {
      mid <- (new_lower + new_upper) / 2
      new_lower <- new_upper <- mid
    }

    refined[[param]] <- .desurv_bo_rebuild_bound_entry(
      entry,
      lower = new_lower,
      upper = new_upper
    )
  }

  refined
}

.desurv_bo_normalise_bound_entry <- function(entry, info_row) {
  if (is.list(entry)) {
    lower <- entry$lower
    upper <- entry$upper
  } else if (is.numeric(entry) && length(entry) == 2L) {
    lower <- entry[1]
    upper <- entry[2]
  } else {
    lower <- info_row$lower
    upper <- info_row$upper
  }
  list(lower = lower, upper = upper)
}

.desurv_bo_rebuild_bound_entry <- function(prev_entry, lower, upper) {
  if (is.list(prev_entry)) {
    prev_entry$lower <- lower
    prev_entry$upper <- upper
    prev_entry
  } else {
    c(lower, upper)
  }
}

.desurv_bo_extract_importance <- function(km_fit, bound_info) {
  params <- bound_info$parameter
  imp <- rep(1, length(params))
  names(imp) <- params

  if (inherits(km_fit, "km")) {
    range_vals <- tryCatch(
      {
        rv <- km_fit@covariance@range.val
        as.numeric(rv)
      },
      error = function(e) NULL
    )
    if (!is.null(range_vals) && length(range_vals) >= length(params)) {
      imp <- 1 / pmax(range_vals[seq_along(params)], .Machine$double.eps)
      imp <- imp / max(imp)
      names(imp) <- params
    }
  }
  imp
}

#' Diagnostics for multi-run Bayesian optimisation refinements
#'
#' @param object Output of [desurv_cv_bayesopt_refine()] or a list of
#'   individual BO runs.
#'
#' @return Invisibly returns a list containing a run-level summary and the
#'   combined history.
#'
#' @export
desurv_bo_refine_diagnostics <- function(object) {
  if (inherits(object, "desurv_cv_bo_refine")) {
    runs <- object$runs
    combined_history <- object$history
    stop_reason <- object$stop_reason
  } else if (is.list(object) && all(vapply(object, inherits, logical(1), "desurv_cv_bo"))) {
    runs <- object
    combined_history <- do.call(rbind, lapply(seq_along(runs), function(i) {
      df <- runs[[i]]$history
      df$run_id <- i
      df
    }))
    stop_reason <- NA_character_
  } else {
    stop("`object` must be a desurv_cv_bo_refine result or list of desurv_cv_bo runs.",
         call. = FALSE)
  }

  n_eval <- vapply(runs, function(run) nrow(run$history), integer(1))
  best <- vapply(runs, function(run) run$best$mean_cindex, numeric(1))
  gain <- c(best[1], diff(best))
  summary_df <- data.frame(
    run_id = seq_along(runs),
    evaluations = n_eval,
    best_mean_cindex = best,
    gain = gain,
    stringsAsFactors = FALSE
  )

  cat("desurv_cv_bayesopt refinement summary\n")
  print(summary_df, row.names = FALSE)
  if (!is.na(stop_reason)) {
    cat("\nStop reason:", stop_reason, "\n")
  }
  cat("\nOverall best mean C-index:",
      sprintf("%.4f", max(best, na.rm = TRUE)), "\n")

  invisible(list(summary = summary_df, history = combined_history))
}

#' Sequentially refine the BO search space for desurv_cv()
#'
#' @inheritParams desurv_cv_bayesopt
#' @param coarse_bounds Initial search ranges supplied to the coarse BO run.
#' @param bo_fixed Named list of fixed hyperparameters (forwarded to
#'   [desurv_cv_bayesopt()]).
#' @param max_refinements Maximum number of refinement rounds after the coarse
#'   run.
#' @param tol_gain Minimum improvement in best mean C-index required to
#'   continue refining.
#' @param plateau Runs without improvement before stopping early.
#' @param top_k Number of top evaluations to define the refined envelope.
#' @param shrink_base Baseline shrinkage factor passed to
#'   [desurv_bo_refine_bounds()].
#' @param importance_gain Additional shrink factor applied using parameter
#'   importances.
#' @param coarse_control Optional named list of arguments overriding the coarse
#'   BO call (e.g., `list(n_init = 10, n_iter = 20)`).
#' @param refine_control Optional named list applied to each refined BO call.
#'
#' @return An object of class `desurv_cv_bo_refine` containing every run,
#'   combined history, the final bounds and the overall best parameters.
#'
#' @export
desurv_cv_bayesopt_refine <- function(
    X,
    y,
    d,
    dataset = NULL,
    samp_keeps = NULL,
    preprocess = TRUE,
    method_trans_train = c("rank", "quant", "none"),
    engine = c("cold", "warmstart"),
    coarse_bounds,
    bo_fixed = list(),
    max_refinements = 3,
    tol_gain = 0.002,
    plateau = 1,
    top_k = 10,
    shrink_base = 0.5,
    importance_gain = 0.3,
    coarse_control = list(),
    refine_control = list(),
    verbose = TRUE,
    ...
) {
  method_trans_train <- match.arg(method_trans_train)
  engine <- match.arg(engine)

  if (!requireNamespace("DiceKriging", quietly = TRUE) ||
      !requireNamespace("lhs", quietly = TRUE)) {
    stop("Packages 'DiceKriging' and 'lhs' are required for desurv_cv_bayesopt_refine().",
         call. = FALSE)
  }

  common_args <- list(
    X = X,
    y = y,
    d = d,
    dataset = dataset,
    samp_keeps = samp_keeps,
    preprocess = preprocess,
    method_trans_train = method_trans_train,
    engine = engine,
    bo_fixed = bo_fixed,
    verbose = verbose,
    ...
  )

  coarse_args <- utils::modifyList(common_args, coarse_control)
  coarse_args$bo_bounds <- coarse_bounds
  coarse_run <- do.call(desurv_cv_bayesopt, coarse_args)

  runs <- list(coarse_run)
  bounds_history <- list(coarse_bounds)
  best_values <- coarse_run$best$mean_cindex
  best_params <- coarse_run$best$params

  no_gain_count <- 0L
  reason <- NA_character_

  for (refine_idx in seq_len(max_refinements)) {
    prev_run <- runs[[length(runs)]]
    prev_bounds <- bounds_history[[length(bounds_history)]]
    imp <- .desurv_bo_extract_importance(prev_run$km_fit, prev_run$bounds)
    refined_bounds <- desurv_bo_refine_bounds(
      prev_bounds,
      history = prev_run$history,
      top_k = top_k,
      importance = imp,
      shrink_base = shrink_base,
      importance_gain = importance_gain
    )

    refine_args <- utils::modifyList(common_args, refine_control)
    refine_args$bo_bounds <- refined_bounds
    refine_args$bo_history <- prev_run
    refined_run <- do.call(desurv_cv_bayesopt, refine_args)

    runs[[length(runs) + 1L]] <- refined_run
    bounds_history[[length(bounds_history) + 1L]] <- refined_bounds

    new_best <- refined_run$best$mean_cindex
    best_values <- c(best_values, new_best)
    best_params <- refined_run$best$params

    gain <- new_best - max(best_values[-length(best_values)])
    if (is.na(gain) || gain < tol_gain) {
      no_gain_count <- no_gain_count + 1L
    } else {
      no_gain_count <- 0L
    }
    if (no_gain_count >= plateau) {
      reason <- sprintf("Stopped after %d refinements with limited improvement.", refine_idx)
      break
    }
  }
  if (is.na(reason)) {
    if (length(runs) - 1L >= max_refinements) {
      reason <- "Reached `max_refinements`."
    } else {
      reason <- "Refinement halted early."
    }
  }

  all_history <- Map(function(run, idx) {
    df <- run$history
    df$run_id <- idx
    df
  }, runs, seq_along(runs))
  combined_history <- do.call(rbind, all_history)
  rownames(combined_history) <- NULL

  overall_best_idx <- which.max(vapply(runs, function(run) run$best$mean_cindex, numeric(1)))
  overall_best <- runs[[overall_best_idx]]$best

  structure(
    list(
      runs = runs,
      bounds = bounds_history,
      history = combined_history,
      overall_best = overall_best,
      stop_reason = reason,
      call = match.call()
    ),
    class = "desurv_cv_bo_refine"
  )
}
