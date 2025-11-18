#' Diagnostics for Bayesian optimisation runs
#'
#' Summarise the evaluation history returned by [desurv_bayesopt()] or
#' [desurv_cv_bayesopt()] and create lightweight diagnostic plots, including a
#' hyperparameter-sensitivity proxy derived from the final Gaussian-process
#' surrogate.
#'
#' @param bo_result Object produced by [desurv_bayesopt()] or
#'   [desurv_cv_bayesopt()].
#' @param include_plots Logical; if `TRUE`, return `ggplot2` objects for
#'   inspection. Plots are skipped when `ggplot2` is not installed.
#'
#' @return A list with components:
#' \itemize{
#'   \item `summary`: tibble-like single-row data frame of aggregate metrics
#'         (evaluation counts, runtime, best score).
#'   \item `progress`: data frame with the cumulative best score per evaluation.
#'   \item `sensitivity`: data frame describing relative importance estimates
#'         derived from the GP length-scales. `NULL` when the surrogate was
#'         unavailable.
#'   \item `plots`: named list of `ggplot` objects (or `NULL` when plotting is
#'         disabled).
#' }
#'
#' @examples
#' \dontrun{
#' bo_fit <- desurv_bayesopt(...)
#' diag <- desurv_bo_diagnostics(bo_fit)
#' diag$plots$progress
#' diag$plots$sensitivity
#' }
#'
#' @export
desurv_bo_diagnostics <- function(bo_result, include_plots = TRUE) {
  if (!inherits(bo_result, c("desurv_bo", "desurv_cv_bo"))) {
    stop("`bo_result` must be the output of desurv_bayesopt() or desurv_cv_bayesopt().",
         call. = FALSE)
  }

  history <- bo_result$history
  if (!is.data.frame(history) || !nrow(history)) {
    stop("`bo_result$history` is not available. Re-run the optimiser with `save_history = TRUE`.",
         call. = FALSE)
  }

  eval_id <- if ("eval_id" %in% names(history)) history$eval_id else seq_len(nrow(history))
  history$eval_id <- eval_id
  history <- history[order(history$eval_id), , drop = FALSE]

  mean_vals <- history$mean_cindex
  status_vec <- if ("status" %in% names(history)) history$status else rep(NA_character_, nrow(history))
  stage_vec <- if ("stage" %in% names(history)) history$stage else rep(NA_character_, nrow(history))
  ok_idx <- status_vec == "ok" & is.finite(mean_vals)
  best_idx <- if (any(ok_idx)) which.max(mean_vals[ok_idx]) else NA_integer_
  best_row <- if (is.na(best_idx)) NULL else history[which(ok_idx)[best_idx], , drop = FALSE]

  cum_best <- suppressWarnings(cummax(ifelse(is.finite(mean_vals), mean_vals, -Inf)))
  cum_best[is.infinite(cum_best)] <- NA_real_
  progress_df <- data.frame(
    eval_id = history$eval_id,
    iteration = if ("iteration" %in% names(history)) history$iteration else NA_integer_,
    mean_cindex = mean_vals,
    best_so_far = cum_best,
    status = status_vec,
    stage = stage_vec,
    stringsAsFactors = FALSE
  )

  total_time <- sum(history$elapsed, na.rm = TRUE)
  summary_row <- data.frame(
    n_evaluations = nrow(history),
    n_success = sum(ok_idx),
    n_errors = nrow(history) - sum(ok_idx),
    total_runtime = total_time,
    best_mean_cindex = if (!is.null(best_row)) best_row$mean_cindex else NA_real_,
    best_eval_id = if (!is.null(best_row)) best_row$eval_id else NA_integer_,
    stringsAsFactors = FALSE
  )

  sensitivity_df <- .desurv_bo_extract_sensitivity(bo_result)

  plot_list <- NULL
  if (isTRUE(include_plots)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' is required for plotting diagnostics. Returning NULL plots.",
              call. = FALSE)
    } else {
      plot_list <- list()
      plot_list$progress <- ggplot2::ggplot(
        progress_df,
        ggplot2::aes(x = eval_id)
      ) +
        ggplot2::geom_line(
          ggplot2::aes(y = best_so_far),
          colour = "#1b7837",
          linewidth = 0.5
        ) +
        ggplot2::geom_point(
          ggplot2::aes(y = mean_cindex, colour = stage, shape = status),
          alpha = 0.7
        ) +
        ggplot2::labs(
          x = "Evaluation",
          y = "Mean C-index",
          title = "Bayesian optimisation progress"
        ) +
        ggplot2::theme_minimal()

      param_cols <- intersect(
        if (!is.null(bo_result$bounds)) bo_result$bounds$parameter else character(),
        names(history)
      )
      if (length(param_cols)) {
        param_df <- data.frame(
          parameter = rep(param_cols, each = nrow(history)),
          value = as.vector(as.matrix(history[, param_cols, drop = FALSE])),
          mean_cindex = rep(history$mean_cindex, times = length(param_cols)),
          stage = rep(stage_vec, times = length(param_cols)),
          status = rep(status_vec, times = length(param_cols)),
          stringsAsFactors = FALSE
        )
        param_df <- param_df[is.finite(param_df$mean_cindex), , drop = FALSE]
        plot_list$parameters <- ggplot2::ggplot(
          param_df,
          ggplot2::aes(x = value, y = mean_cindex, colour = stage)
        ) +
          ggplot2::geom_point(alpha = 0.5, size = 1) +
          ggplot2::facet_wrap(~parameter, scales = "free_x") +
          ggplot2::labs(
            x = "Hyperparameter value",
            y = "Mean C-index",
            title = "Observed parameter effects"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position = "bottom")
      }

      if (!is.null(sensitivity_df)) {
        plot_list$sensitivity <- ggplot2::ggplot(
          sensitivity_df,
          ggplot2::aes(
            x = stats::reorder(parameter, importance),
            y = importance
          )
        ) +
          ggplot2::geom_col(fill = "#2166ac") +
          ggplot2::coord_flip() +
          ggplot2::labs(
            x = "Hyperparameter",
            y = "Relative importance",
            title = "GP-based sensitivity"
          ) +
          ggplot2::theme_minimal()
      }
    }
  }

  structure(
    list(
      summary = summary_row,
      progress = progress_df,
      sensitivity = sensitivity_df,
      plots = plot_list,
      call = match.call()
    ),
    class = "desurv_bo_diagnostics"
  )
}

.desurv_bo_extract_sensitivity <- function(bo_result) {
  km_fit <- bo_result$km_fit
  if (is.null(km_fit)) {
    return(NULL)
  }
  lengthscale <- tryCatch(
    km_fit@covariance@range.val,
    error = function(e) NULL
  )
  if (is.null(lengthscale)) {
    lengthscale <- tryCatch(km_fit@covariance@range.par, error = function(e) NULL)
  }
  if (is.null(lengthscale)) {
    return(NULL)
  }

  params <- if (!is.null(bo_result$bounds)) bo_result$bounds$parameter else names(bo_result$best$params)
  params <- params[seq_along(lengthscale)]
  lengthscale <- lengthscale[seq_along(params)]
  inv_scale <- 1 / pmax(lengthscale, .Machine$double.eps)
  rel_imp <- inv_scale / sum(inv_scale)
  data.frame(
    parameter = params,
    lengthscale = lengthscale,
    importance = rel_imp,
    stringsAsFactors = FALSE
  )
}
