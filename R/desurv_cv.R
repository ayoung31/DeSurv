#' Cross-validated DeSurv model selection and refitting
#'
#' @description
#' Perform cross-validation over a grid of hyperparameters for the DeSurv
#' model, select the combination that maximizes mean validation C-index
#' (or applies a 1-SE rule), and refit a final model on the full dataset
#' using the chosen hyperparameters.
#'
#' @details
#' This function is a high-level wrapper around the internal workhorse
#' \code{.desurv_cv_grid_warmstart()}, which:
#' \itemize{
#'   \item builds a grid of hyperparameters
#'         \code{k, lambda, nu, lambdaW, lambdaH},
#'   \item performs K-fold cross-validation for each combination, using
#'         warmstarts along \code{alpha_grid} within each training fold via
#'         \code{\link{desurv_alpha_warmstart}},
#'   \item computes validation C-indices for each
#'         \code{(k, lambda, nu, lambdaW, lambdaH, alpha, fold, init_id)},
#'   \item summarizes results in two stages:
#'         \enumerate{
#'           \item For each fold, averages the C-index across random
#'                 initialization paths (inits), yielding a fold-level mean
#'                 C-index.
#'           \item Across folds, computes the overall mean and standard error
#'                 (SE) of these fold-level means. The SE is calculated as the
#'                 standard deviation of the fold-level means divided by the
#'                 square root of the number of non-missing folds.
#'         }
#' }
#'
#' The present function uses these precomputed summaries directly from
#' \code{.desurv_cv_grid_warmstart()}:
#' \itemize{
#'   \item \code{summary_fold}: contains the mean C-index across inits for
#'         each fold, \code{alpha}, and hyperparameter combination.
#'   \item \code{summary}: contains the mean and SE of C-index across folds
#'         (after averaging across inits within each fold) for each
#'         \code{(k, lambda, nu, lambdaW, lambdaH, alpha)} combination.
#' }
#'
#' The \code{rule} argument controls the selection criterion:
#'
#' \describe{
#'   \item{\code{"max"}}{
#'     Selects the combination with the largest mean C-index.
#'   }
#'   \item{\code{"1se"}}{
#'     Let \eqn{\hat{c}^*} denote the largest mean C-index and
#'     \eqn{\mathrm{SE}^*} its standard error. Define a threshold
#'     \eqn{\hat{c}^* - \mathrm{SE}^*}. Among all combinations with mean
#'     C-index at least this threshold, select the "simplest" one using the
#'     deterministic tie-breaking order below. If the SE at the best
#'     combination is unavailable (e.g., fewer than 2 folds), the rule
#'     falls back to \code{"max"}.
#'   }
#' }
#'
#' Ties between equivalent mean C-indices are broken deterministically in
#' the following order:
#' \itemize{
#'   \item smallest \code{k},
#'   \item then smallest \code{alpha},
#'   \item then largest \code{lambda},
#'   \item then smallest \code{nu},
#'   \item then smallest \code{lambdaW},
#'   \item then smallest \code{lambdaH}.
#' }
#'
#' After selecting the optimal hyperparameters, a final model is refit on the
#' full dataset by calling \code{\link{desurv_fit}} with the chosen settings.
#'
#' @param X Either a numeric matrix (\eqn{p \times n}) of nonnegative features
#'   (e.g., genes by samples), or a \code{"desurv_data"} object created by
#'   \code{\link{desurv_data}}.
#' @param y Numeric vector of survival times of length \eqn{n}. Used only if
#'   \code{X} is a matrix or data frame; ignored if \code{X} is a
#'   \code{"desurv_data"} object.
#' @param d Numeric (0/1) vector of event indicators of length \eqn{n}.
#'   Used only if \code{X} is a matrix or data frame; ignored if \code{X}
#'   is a \code{"desurv_data"} object.
#'
#' @param k_grid Numeric or integer vector of candidate ranks (number of latent
#'   factors) to consider.
#' @param alpha_grid Numeric vector of supervision parameter values in
#'   \eqn{[0, 1]} to explore during CV.
#' @param lambda_grid Numeric vector of candidate global penalty parameters.
#' @param nu_grid Numeric vector of candidate elastic-net mixing parameters
#'   in \eqn{[0, 1]}.
#' @param lambdaW_grid Numeric vector of candidate penalties for the
#'   \code{W} factor.
#' @param lambdaH_grid Numeric vector of candidate penalties for the
#'   \code{H} factor.
#'
#' @param n_starts Integer; number of random initialization paths per
#'   hyperparameter combination and fold. Each path is warmstarted along
#'   \code{alpha_grid} via \code{\link{desurv_alpha_warmstart}}.
#' @param nfolds Integer; number of CV folds.
#' @param folds Optional integer vector of length \code{n} with values in
#'   \code{1:nfolds} specifying the fold assignment of each sample. If
#'   \code{NULL}, stratified folds are generated automatically so that
#'   events and non-events are approximately balanced.
#' @param seed Integer or \code{NULL}; base seed used for fold generation
#'   and for randomness inside the warmstart paths.
#'
#' @param tol Numeric convergence tolerance for the full optimizations inside
#'   \code{\link{desurv_fit}} during CV and in the final refit.
#' @param maxit Integer maximum number of iterations for each full
#'   optimization inside \code{\link{desurv_fit}} during CV and in the final
#'   refit.
#'
#' @param rule Character string specifying the model selection rule; one of
#'   \code{"max"} (default) or \code{"1se"}. See Details.
#'
#' @param parallel_grid Logical; if \code{TRUE}, the CV computations are
#'   parallelized over the (hyperparameter, fold) jobs using
#'   \code{parallel::mclapply()} on non-Windows platforms.
#' @param ncores_grid Integer or \code{NULL}; number of cores to use when
#'   \code{parallel_grid = TRUE}. If \code{NULL}, defaults to
#'   \code{parallel::detectCores()} capped at the number of such jobs.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed.
#'
#' @return
#' An object of class \code{"desurv_cv"} with components:
#' \itemize{
#'   \item \code{fit}: the final \code{"desurv_fit"} object refit on the full
#'         data using the chosen hyperparameter combination.
#'   \item \code{best}: a list describing the selected hyperparameters, with
#'         elements \code{k}, \code{lambda}, \code{nu}, \code{lambdaW},
#'         \code{lambdaH}, \code{alpha}, \code{mean_cindex}, and
#'         \code{se_cindex}.
#'   \item \code{rule}: the selection rule used (\code{"max"} or \code{"1se"}).
#'   \item \code{cv_results}: data frame of raw fold-level CV results with
#'         columns \code{fold}, \code{k}, \code{lambda}, \code{nu},
#'         \code{lambdaW}, \code{lambdaH}, \code{init_id}, \code{alpha},
#'         and \code{cindex}.
#'   \item \code{summary_fold}: data frame containing mean C-index across
#'         inits for each fold, alpha, and hyperparameter combination.
#'   \item \code{summary}: data frame containing the mean and SE of the
#'         C-index across folds (after averaging across inits within each
#'         fold) for each \code{(hyperparameters, alpha)}.
#'   \item \code{alpha_grid}: the alpha grid used.
#'   \item \code{hyper_grid}: the hyperparameter grid (excluding alpha).
#'   \item \code{folds}: the fold assignments used.
#'   \item \code{call}: the original function call.
#' }
#'
#' @seealso
#' \code{\link{desurv_fit}}, \code{\link{desurv_alpha_warmstart}},
#' \code{\link{.desurv_cv_grid_warmstart}}
#'
#' @export
desurv_cv <- function(
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
    rule          = c("max", "1se"),
    parallel_grid = FALSE,
    ncores_grid   = NULL,
    verbose       = TRUE
) {
  call <- match.call()
  rule <- match.arg(rule)

  # Run internal CV engine (now returns: cv_results, summary_fold, summary)
  cv_obj <- .desurv_cv_grid_warmstart(
    X            = X,
    y            = y,
    d            = d,
    k_grid       = k_grid,
    alpha_grid   = alpha_grid,
    lambda_grid  = lambda_grid,
    nu_grid      = nu_grid,
    lambdaW_grid = lambdaW_grid,
    lambdaH_grid = lambdaH_grid,
    n_starts      = n_starts,
    nfolds        = nfolds,
    folds         = folds,
    seed          = seed,
    tol           = tol,
    maxit         = maxit,
    parallel_grid = parallel_grid,
    ncores_grid   = ncores_grid,
    verbose       = verbose
  )

  # Basic sanity check
  if (is.null(cv_obj$summary) || nrow(cv_obj$summary) == 0L) {
    stop("No CV summary available to select a best model.")
  }

  # Use the pre-computed hyper+alpha summaries (mean + SE across folds)
  summ2 <- cv_obj$summary

  mc <- summ2$mean_cindex
  se <- summ2$se_cindex

  # tie-breaking order: k (small), alpha (small), lambda (large), etc.
  tie_order <- function(df) {
    with(
      df,
      order(
        k,           # smallest k
        alpha,       # then smallest alpha
        -lambda,     # then largest lambda
        nu,          # then smallest nu
        lambdaW,     # then smallest lambdaW
        lambdaH      # then smallest lambdaH
      )
    )
  }

  # Choose best hyperparameters
  if (rule == "max") {
    max_c <- max(mc, na.rm = TRUE)
    cand  <- summ2[mc == max_c, , drop = FALSE]
    best_row <- cand[tie_order(cand)[1L], , drop = FALSE]

  } else { # "1se"
    idx_best <- which.max(mc)
    best0    <- summ2[idx_best, , drop = FALSE]
    se0      <- best0$se_cindex

    if (is.na(se0)) {
      # fallback to max rule if SE not available
      max_c <- max(mc, na.rm = TRUE)
      cand  <- summ2[mc == max_c, , drop = FALSE]
      best_row <- cand[tie_order(cand)[1L], , drop = FALSE]
    } else {
      thr  <- best0$mean_cindex - se0
      cand <- summ2[mc >= thr, , drop = FALSE]
      if (nrow(cand) == 0L) {
        max_c <- max(mc, na.rm = TRUE)
        cand  <- summ2[mc == max_c, , drop = FALSE]
      }
      best_row <- cand[tie_order(cand)[1L], , drop = FALSE]
    }
  }

  if (verbose) {
    msg_rule <- if (rule == "max") "max-mean rule" else "1-SE rule"
    message(sprintf(
      "Selected (rule = %s): k = %d, lambda = %.3g, nu = %.3g, lambdaW = %.3g, lambdaH = %.3g, alpha = %.3g (mean CV c-index = %.4f).",
      msg_rule,
      best_row$k, best_row$lambda, best_row$nu,
      best_row$lambdaW, best_row$lambdaH,
      best_row$alpha, best_row$mean_cindex
    ))
  }

  # Refit final model on full data with best hyperparameters
  if (inherits(X, "desurv_data")) {
    fit_best <- desurv_fit(
      X       = X,
      alpha   = best_row$alpha,
      lambda  = best_row$lambda,
      nu      = best_row$nu,
      lambdaW = best_row$lambdaW,
      lambdaH = best_row$lambdaH,
      tol     = tol,
      maxit   = maxit,
      verbose = verbose
    )
  } else {
    fit_best <- desurv_fit(
      X       = X,
      y       = y,
      d       = d,
      k       = best_row$k,
      alpha   = best_row$alpha,
      lambda  = best_row$lambda,
      nu      = best_row$nu,
      lambdaW = best_row$lambdaW,
      lambdaH = best_row$lambdaH,
      tol     = tol,
      maxit   = maxit,
      verbose = verbose
    )
  }

  best_list <- list(
    k           = best_row$k,
    lambda      = best_row$lambda,
    nu          = best_row$nu,
    lambdaW     = best_row$lambdaW,
    lambdaH     = best_row$lambdaH,
    alpha       = best_row$alpha,
    mean_cindex = best_row$mean_cindex,
    se_cindex   = best_row$se_cindex
  )

  out <- list(
    fit          = fit_best,
    best         = best_list,
    rule         = rule,
    cv_results   = cv_obj$cv_results,   # raw per fold, per init, per alpha, per hyper
    summary_fold = cv_obj$summary_fold, # mean across inits, per fold+alpha+hyper
    summary      = cv_obj$summary,      # mean across folds + SE, per alpha+hyper
    alpha_grid   = cv_obj$alpha_grid,
    hyper_grid   = cv_obj$hyper_grid,
    folds        = cv_obj$folds,
    call         = call
  )
  class(out) <- "desurv_cv"
  out
}

