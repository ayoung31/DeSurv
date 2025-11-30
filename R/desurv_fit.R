#' Fit a DeSurv model
#'
#' @description
#' Fit a survival-driven nonnegative matrix factorization (DeSurv) model that
#' couples an NMF decomposition of \code{X} with a Cox proportional hazards
#' model on the latent factors. The function supports both automatic
#' multi-start initialization and user-specified initial values for
#' \code{W}, \code{H}, and \code{beta}.
#'
#' @details
#' The core model decomposes a nonnegative data matrix \eqn{X \in R^{p \times n}}
#' into nonnegative factors \eqn{W \in R^{p \times k}} and \eqn{H \in R^{k \times n}},
#' and fits a Cox model with linear predictor
#' \deqn{\eta = (X^\top W)\,\beta,}
#' where \eqn{\beta \in R^k}. The parameters are:
#' \itemize{
#'   \item \code{alpha}: supervision parameter in \eqn{[0, 1]} controlling the
#'         strength of the survival-driven component of the objective.
#'   \item \code{nu}: elastic-net mixing parameter in \eqn{[0, 1]} for the
#'         penalty on \eqn{\beta}, interpolating between ridge-like and
#'         lasso-like behaviour.
#' }
#'
#' \strong{Initialization strategy}
#'
#' \itemize{
#'   \item If \strong{no} user-specified \code{W0}, \code{H0}, \code{beta0} are
#'         supplied, the function first runs \code{ninit} short optimizations
#'         via an internal \code{init()} routine, each with at most
#'         \code{imaxit} iterations and tolerance \code{tol_init}. These short
#'         runs can be executed in parallel if \code{parallel_init = TRUE}.
#'         The best partial fit is selected by concordance index (C-index),
#'         and a single full optimization is then performed from that
#'         initialization using \code{maxit} iterations and tolerance
#'         \code{tol}.
#'   \item If the user supplies \code{W0} (with optional \code{H0} and
#'         \code{beta0}), these are used directly as the starting point for a
#'         single full optimization with \code{maxit} iterations and tolerance
#'         \code{tol}; in this case, \code{tol_init}, \code{imaxit},
#'         \code{ninit}, \code{parallel_init}, and \code{ncores_init} are
#'         \strong{ignored}. When \code{H0} is omitted it is derived from
#'         \code{W0} and the training matrix via non-negative least squares,
#'         and \code{beta0} defaults to a zero vector.
#' }
#'
#' The function accepts either raw inputs \code{X, y, d, k} or a pre-constructed
#' \code{"desurv_data"} object created by \code{\link{desurv_data}} as the first
#' argument.
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
#' @param k Positive integer rank (number of latent factors). Used only if
#'   \code{X} is a matrix or data frame; ignored if \code{X} is a
#'   \code{"desurv_data"} object (where \code{k} is already stored).
#' @param alpha Numeric scalar in \eqn{[0, 1]} controlling the supervision
#'   strength (weight of the survival component in the objective).
#' @param lambda Numeric non-negative scalar; global penalty parameter.
#' @param nu Numeric scalar in \eqn{[0, 1]} giving the elastic-net mixing
#'   parameter for \eqn{\beta}.
#' @param lambdaW Numeric non-negative scalar; penalty parameter applied to
#'   the \code{W} factor.
#' @param lambdaH Numeric non-negative scalar; penalty parameter applied to
#'   the \code{H} factor.
#' @param W0 Optional numeric matrix of dimension \eqn{p \times k} giving a
#'   user-specified initial value for \code{W}. Providing \code{W0} skips the
#'   internal multi-start phase.
#' @param H0 Optional numeric matrix of dimension \eqn{k \times n} giving a
#'   user-specified initial value for \code{H}. When omitted (but \code{W0} is
#'   supplied) \code{H0} is computed internally via non-negative least squares.
#' @param beta0 Optional numeric vector of length \code{k} giving a
#'   user-specified initial value for the Cox coefficients \code{beta}.
#'   Defaults to a zero vector when omitted in combination with a custom
#'   \code{W0}.
#' @param seed Integer or \code{NULL}; random seed used during initialization
#'   (for the short runs when \code{W0, H0, beta0} are not provided).
#' @param tol Numeric convergence tolerance for the \strong{full} optimization
#'   run. This is passed to the underlying optimizer for the final fit.
#' @param tol_init Numeric convergence tolerance used during the \strong{short
#'   initialization runs} inside \code{init()}. Ignored if user-specified
#'   \code{W0, H0, beta0} are provided.
#' @param maxit Integer giving the maximum number of iterations for the full
#'   optimization starting from the chosen initialization.
#' @param imaxit Integer giving the maximum number of iterations for each
#'   short initialization run in \code{init()}. Ignored if user-specified
#'   \code{W0, H0, beta0} are provided.
#' @param ninit Integer number of random initializations to try in
#'   \code{init()}. Ignored if user-specified \code{W0, H0, beta0} are
#'   provided.
#' @param parallel_init Logical; if \code{TRUE}, the short initialization runs
#'   in \code{init()} are executed in parallel using
#'   \code{parallel::mclapply()} on non-Windows platforms. Ignored if
#'   user-specified \code{W0, H0, beta0} are provided.
#' @param ncores_init Optional integer specifying the number of cores to use
#'   when \code{parallel_init = TRUE}. If \code{NULL}, defaults to
#'   \code{parallel::detectCores()} capped at \code{ninit}. Ignored if
#'   user-specified \code{W0, H0, beta0} are provided.
#' @param verbose Logical; if \code{TRUE}, progress messages and warnings are
#'   printed during initialization and optimization.
#' @param preprocess_info Optional list containing preprocessing metadata
#'   (gene order, transformation method, normalization targets) that should
#'   be stored on the fitted object for use during prediction.
#'
#' @return
#' An object of class \code{"desurv_fit"} with components including:
#' \itemize{
#'   \item \code{W}: estimated \eqn{p \times k} basis matrix.
#'   \item \code{H}: estimated \eqn{k \times n} coefficient matrix.
#'   \item \code{beta}: estimated Cox coefficients of length \code{k}.
#'   \item \code{cindex}: C-index of the final fit on the training data.
#'   \item \code{cindex_init}: C-index of the best short initialization
#'         (if multi-start was used), or \code{NA} when user-specified
#'         initialization is provided.
#'   \item \code{convergence}: logical flag indicating whether the final
#'         optimization satisfied the stopping criterion before exhausting
#'         \code{maxit}.
#'   \item \code{data}: the underlying \code{"desurv_data"} object used for
#'         fitting.
#'   \item \code{hyper}: list of hyperparameters actually used after
#'         validation and possible clipping.
#'   \item \code{preprocess}: optional preprocessing metadata captured from
#'         the training pipeline.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' p <- 50
#' n <- 100
#' k <- 3
#'
#' X <- matrix(rexp(p * n), nrow = p, ncol = n)
#' y <- rexp(n)
#' d <- rbinom(n, 1, 0.7)
#'
#' # Basic usage with raw inputs and automatic multi-start
#' fit <- desurv_fit(
#'   X, y, d, k = k,
#'   alpha   = 0.5,
#'   lambda  = 1,
#'   nu      = 0.5,
#'   lambdaW = 0.1,
#'   lambdaH = 0.1,
#'   tol     = 1e-6,
#'   tol_init = 1e-4,
#'   maxit   = 500,
#'   imaxit  = 20,
#'   ninit   = 5,
#'   parallel_init = TRUE,
#'   verbose = TRUE
#' )
#'
#' # Using a pre-constructed desurv_data object
#' dd  <- desurv_data(X, y, d, k)
#' fit2 <- desurv_fit(
#'   dd,
#'   alpha   = 0.5,
#'   lambda  = 1,
#'   nu      = 0.5,
#'   lambdaW = 0.1,
#'   lambdaH = 0.1
#' )
#'
#' # Using user-specified initialization (skips multi-start)
#' W0 <- matrix(runif(p * k, 0, max(X)), nrow = p)
#' H0 <- matrix(runif(n * k, 0, max(X)), nrow = k)
#' beta0 <- rep(0, k)
#'
#' fit3 <- desurv_fit(
#'   X, y, d, k = k,
#'   alpha   = 0.5,
#'   lambda  = 1,
#'   nu      = 0.5,
#'   lambdaW = 0.1,
#'   lambdaH = 0.1,
#'   W0 = W0, H0 = H0, beta0 = beta0,
#'   maxit   = 500,
#'   tol     = 1e-6
#' )
#' }
#'
#' @export
desurv_fit <- function(
    X, y = NULL, d = NULL, k = NULL,
    alpha, lambda, nu, lambdaW, lambdaH,
    W0            = NULL,
    H0            = NULL,
    beta0         = NULL,
    seed          = 123,
    tol           = 1e-6,
    tol_init      = 1e-4,
    maxit         = 1000,
    imaxit        = 30,
    ninit         = 20,
    parallel_init = FALSE,
    ncores_init   = NULL,
    verbose       = TRUE,
    preprocess_info = NULL
) {
  # Build or reuse desurv_data object
  data <- .as_desurv_data(X, y = y, d = d, k = k)

  gene_names <- rownames(data$X)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(data$X)))
    rownames(data$X) <- gene_names
  }
  attach_gene_names <- function(mat) {
    if (is.null(mat)) {
      return(mat)
    }
    mat <- as.matrix(mat)
    if (nrow(mat) != length(gene_names)) {
      stop("Internal error: gene dimension mismatch in initialization.", call. = FALSE)
    }
    if (is.null(rownames(mat))) {
      rownames(mat) <- gene_names
    }
    mat
  }

  # Validate hyperparameters (alpha, lambda, nu, lambdaW, lambdaH, tol, maxit, verbose)
  hp <- .validate_desurv_hyperparams(
    alpha, lambda, nu, lambdaW, lambdaH,
    tol, maxit, verbose
  )

  # light check for imaxit
  imaxit_num <- as.numeric(imaxit)
  if (length(imaxit_num) != 1L || is.na(imaxit_num) || imaxit_num <= 0) {
    warning("imaxit not positive; setting imaxit = 30.")
    imaxit <- 30L
  } else {
    imaxit <- as.integer(imaxit_num)
  }

  # Decide whether we're using a user-specified initialization
  use_custom_init <- !is.null(W0) || !is.null(H0) || !is.null(beta0)

  if (use_custom_init) {
    init_par <- .validate_desurv_custom_init(data$X, data$k, W0, H0, beta0)
    best_W   <- attach_gene_names(init_par$W0)
    best_H   <- init_par$H0
    best_beta <- init_par$beta0
    cindex_init <- NA_real_
  } else {
    # Multi-start short initializations to pick a good starting point
    init_res <- init(
      X       = data$X,
      y       = data$y,
      delta   = data$d,
      k       = data$k,
      alpha   = hp$alpha,
      lambda  = hp$lambda,
      nu      = hp$nu,
      lambdaW = hp$lambdaW,
      lambdaH = hp$lambdaH,
      seed    = seed,
      tol_init = tol_init,
      imaxit  = imaxit,
      verbose = verbose,
      ninit   = ninit,
      parallel = parallel_init,
      ncores   = ncores_init
    )

    best_W      <- init_res$W0
    best_H      <- init_res$H0
    best_beta   <- init_res$beta0
    cindex_init <- init_res$cindex_init

    best_W <- attach_gene_names(best_W)
  }

  # Optional: reseed before full run (only for auto-init path)
  if (!is.null(seed) && !use_custom_init) {
    set.seed(as.integer(seed))
  }

  # Full optimization from chosen initialization
  fit_full <- .run_optimize_loss(
    data$X, data$y, data$d,
    best_W, best_H, best_beta,
    hp$alpha, hp$lambda, hp$nu,
    hp$lambdaW, hp$lambdaH,
    hp$tol, hp$maxit, verbose
  )

  fit_full$W <- attach_gene_names(fit_full$W)

  Z_full  <- t(data$X) %*% fit_full$W
  lp_full <- drop(Z_full %*% fit_full$beta)
  cindex_full <- cvwrapr::getCindex(lp_full, survival::Surv(data$y, data$d))

  structure(
    list(
      W           = fit_full$W,
      H           = fit_full$H,
      beta        = fit_full$beta,
      cindex      = cindex_full,
      cindex_init = cindex_init,
      convergence = isTRUE(fit_full$convergence),
      data        = data,
      hyper       = hp,
      preprocess  = preprocess_info
    ),
    class = "desurv_fit"
  )
}
