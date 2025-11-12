#' Cross-validated DeSurv model selection and refitting
#'
#' @description
#' Perform cross-validation over a grid of hyperparameters for the DeSurv
#' model, select the combination that maximizes mean validation C-index
#' (or applies a 1-SE rule), and refit a final model on the full dataset
#' using the chosen hyperparameters.
#'
#' @details
#' This function is a high-level wrapper around the internal CV helpers
#' \code{.desurv_cv_grid_warmstart()} (default) and
#' \code{.desurv_cv_grid()} (when \code{engine = "cold"}), which:
#' \itemize{
#'   \item builds a grid of hyperparameters
#'         \code{k, lambda, nu, lambdaW, lambdaH},
#'   \item performs K-fold cross-validation for each combination; with
#'         \code{engine = "warmstart"} the fits are warmstarted along
#'         \code{alpha_grid} inside each fold via
#'         \code{\link{desurv_alpha_warmstart}}, whereas the \code{"cold"}
#'         engine refits each \code{alpha} independently,
#'   \item computes validation C-indices for each
#'         \code{(k, lambda, nu, lambdaW, lambdaH, alpha, fold)} combination.
#'         With the warmstart engine we retain multiple initialization paths
#'         (tracked via \code{init_id}); the cold engine forwards
#'         \code{n_starts} to \code{desurv_fit()} as \code{ninit} and stores
#'         a single score per combination.
#'   \item summarizes results in two stages:
#'         \enumerate{
#'           \item For each fold, averages the C-index across random
#'                 initialization paths when they exist (warmstart engine);
#'                 otherwise the fold-level mean is the single available score.
#'           \item Across folds, computes the overall mean and standard error
#'                 (SE) of these fold-level means. The SE is calculated as the
#'                 standard deviation of the fold-level means divided by the
#'                 square root of the number of non-missing folds.
#'         }
#' }
#'
#' Optionally, users can request that the raw expression matrix be filtered
#' and transformed via \code{\link{preprocess_data}} before CV by setting
#' \code{preprocess = TRUE} and supplying the associated arguments
#' (\code{dataset}, \code{samp_keeps}, \code{ngene}, \code{genes},
#' \code{method_trans_train}). This subsumes the behavior previously exposed
#' via \code{\link{desurv_cv_preprocess}}.
#'
#' The present function uses these precomputed summaries directly from the
#' selected engine:
#' \itemize{
#'   \item \code{summary_fold}: contains the mean C-index for
#'         each fold, \code{alpha}, and hyperparameter combination (after
#'         averaging across initialization paths when applicable).
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
#' @param dataset Optional vector identifying the dataset or cohort for each
#'   sample. When provided (or embedded within a \code{"desurv_data"}
#'   object), automatically generated folds are stratified jointly by
#'   event indicator and dataset.
#' @param preprocess Logical; if \code{TRUE}, run
#'   \code{\link{preprocess_data}} on the raw inputs before cross-validation.
#'   Requires that \code{X} be a numeric matrix/data frame (not a
#'   \code{"desurv_data"} object) and that \code{dataset} be supplied.
#' @param samp_keeps Optional logical/integer/character index of samples to
#'   keep prior to preprocessing (only used when \code{preprocess = TRUE}).
#' @param ngene Integer giving the number of genes to retain per dataset when
#'   \code{preprocess = TRUE} and \code{genes} is \code{NULL}.
#' @param genes Optional character vector of genes to keep exactly when
#'   \code{preprocess = TRUE}. Overrides \code{ngene} selection.
#' @param method_trans_train Character; one of \code{"rank"},
#'   \code{"quant"}, or \code{"none"}, forwarded to
#'   \code{\link{preprocess_data}} when \code{preprocess = TRUE}.
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
#'   hyperparameter combination and fold. When \code{engine = "warmstart"},
#'   each path is warmstarted along \code{alpha_grid} via
#'   \code{\link{desurv_alpha_warmstart}}. The \code{"cold"} engine now
#'   forwards \code{n_starts} directly to \code{\link{desurv_fit}} as
#'   \code{ninit}, so each \code{alpha} inside a fold is refit once using
#'   that many internal starts.
#' @param nfolds Integer; number of CV folds.
#' @param folds Optional integer vector of length \code{n} with values in
#'   \code{1:nfolds} specifying the fold assignment of each sample. If
#'   \code{NULL}, stratified folds are generated automatically so that
#'   events and non-events are approximately balanced.
#' @param seed Integer or \code{NULL}; base seed used for fold generation
#'   and for randomness inside the initialization paths.
#'
#' @param tol Numeric convergence tolerance for the full optimizations inside
#'   \code{\link{desurv_fit}} during CV and in the final refit.
#' @param maxit Integer maximum number of iterations for each full
#'   optimization inside \code{\link{desurv_fit}} during CV and in the final
#'   refit.
#' @param engine Character string specifying which internal CV engine to
#'   use. The \code{"warmstart"} engine reuses alpha paths via
#'   \code{\link{desurv_alpha_warmstart}}; \code{"cold"} calls
#'   \code{\link{.desurv_cv_grid}} to fit each alpha independently.
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
#' @param cv_only Logical; if \code{TRUE}, skip the final refit and return
#'   the raw cross-validation object produced by the selected engine.
#' 
#' @return
#' An object of class \code{"desurv_cv"} with components:
#' \itemize{
#'   \item \code{fit}: the final \code{"desurv_fit"} object refit on the full
#'         data using the chosen hyperparameter combination (omitted when
#'         \code{cv_only = TRUE}).
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
#'   \item \code{preprocess}: metadata describing any preprocessing that was
#'         applied (only present when \code{preprocess = TRUE}).
#'   \item \code{call}: the original function call.
#' }
#'
#' @seealso
#' \code{\link{desurv_fit}}, \code{\link{desurv_alpha_warmstart}},
#' \code{\link{preprocess_data}}, \code{\link{desurv_cv_preprocess}},
#' \code{\link{.desurv_cv_grid_warmstart}}, \code{\link{.desurv_cv_grid}}
#'
#' @export
desurv_cv <- function(
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
    engine        = c("warmstart", "cold"),
    rule          = c("max", "1se"),
    parallel_grid = FALSE,
    ncores_grid   = NULL,
    verbose       = TRUE,
    preprocess    = FALSE,
    samp_keeps    = NULL,
    ngene         = 1000,
    genes         = NULL,
    method_trans_train = c("rank", "quant", "none"),
    cv_only       = FALSE
) {
  call <- match.call()
  rule <- match.arg(rule)
  engine <- match.arg(engine)
  method_trans_train <- match.arg(method_trans_train)

  prepared <- .desurv_cv_prepare_inputs(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    preprocess         = preprocess,
    samp_keeps         = samp_keeps,
    ngene              = ngene,
    genes              = genes,
    method_trans_train = method_trans_train,
    verbose            = verbose
  )

  X <- prepared$X
  y <- prepared$y
  d <- prepared$d
  dataset <- prepared$dataset
  preprocess_info <- prepared$preprocess_info

  cv_fun <- switch(
    engine,
    warmstart = .desurv_cv_grid_warmstart,
    cold      = .desurv_cv_grid
  )

  # Run internal CV engine (now returns: cv_results, summary_fold, summary)
  cv_obj <- cv_fun(
    X            = X,
    y            = y,
    d            = d,
    dataset      = dataset,
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

  if (cv_only) {
    if (!is.null(preprocess_info)) {
      cv_obj$preprocess <- preprocess_info
    }
    cv_obj$rule   <- rule
    cv_obj$engine <- engine
    cv_obj$call   <- call
    cv_obj$fit    <- NULL
    class(cv_obj) <- "desurv_cv"
    return(cv_obj)
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
  if (!is.null(preprocess_info)) {
    out$preprocess <- preprocess_info
  }
  class(out) <- "desurv_cv"
  out
}

.desurv_cv_prepare_inputs <- function(
    X,
    y,
    d,
    dataset,
    preprocess,
    samp_keeps,
    ngene,
    genes,
    method_trans_train,
    verbose
) {
  if (!preprocess) {
    return(list(
      X = X,
      y = y,
      d = d,
      dataset = dataset,
      preprocess_info = NULL
    ))
  }

  if (inherits(X, "desurv_data")) {
    stop("preprocess = TRUE is not supported when `X` is already a 'desurv_data' object.")
  }
  if (is.null(dataset)) {
    stop("`dataset` must be supplied when preprocess = TRUE.")
  }

  if (isTRUE(verbose)) {
    message(paste0("Filtering to top ", ngene, " highly expressed and variable genes per dataset"))
  }

  data_proc <- preprocess_data(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    samp_keeps         = samp_keeps,
    ngene              = ngene,
    genes              = genes,
    method_trans_train = method_trans_train,
    verbose            = verbose
  )

  required_cols <- c("time", "event")
  if (!all(required_cols %in% colnames(data_proc$sampInfo))) {
    stop("Expected columns 'time' and 'event' in data$sampInfo after preprocessing.")
  }

  X_proc <- data_proc$ex
  y_proc <- data_proc$sampInfo$time
  d_proc <- data_proc$sampInfo$event
  dataset_proc <- data_proc$sampInfo$dataset

  genes_proc <- if (!is.null(data_proc$featInfo)) {
    data_proc$featInfo
  } else if (!is.null(rownames(X_proc))) {
    rownames(X_proc)
  } else {
    genes
  }

  preprocess_info <- list(
    ngene              = ngene,
    genes              = genes_proc,
    method_trans_train = method_trans_train,
    samp_keeps         = data_proc$samp_keeps
  )

  list(
    X = X_proc,
    y = y_proc,
    d = d_proc,
    dataset = dataset_proc,
    preprocess_info = preprocess_info
  )
}
