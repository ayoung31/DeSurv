#' @keywords internal
#' @title Stratified fold assignment by event indicator
#' @description
#' Create stratified K-fold assignments so that events (d == 1) and
#' datasets are approximately balanced across folds.
.desurv_make_folds_stratified <- function(d, dataset = NULL, nfolds, seed = NULL) {
  d <- as.numeric(d)
  if (any(is.na(d)) || !all(d %in% c(0, 1))) {
    stop("d must be a 0/1 vector without NA for stratified folds.")
  }

  nfolds <- as.integer(nfolds)
  if (length(nfolds) != 1L || is.na(nfolds) || nfolds <= 1L) {
    stop("nfolds must be an integer >= 2.")
  }

  n <- length(d)
  if (is.null(dataset)) {
    dataset <- rep("_all", n)
  } else {
    if (length(dataset) != n) {
      stop("`dataset` must have length equal to length(d).")
    }
    if (any(is.na(dataset))) {
      stop("`dataset` cannot contain NA.")
    }
  }
  dataset <- as.character(dataset)

  if (!is.null(seed)) set.seed(as.integer(seed))

  strata <- interaction(d, dataset, drop = TRUE)
  folds <- integer(n)
  levels_strata <- levels(strata)

  for (lvl in levels_strata) {
    idx <- which(strata == lvl)
    if (length(idx) == 0L) next
    idx <- sample(idx, length(idx))
    folds_level <- rep(seq_len(nfolds), length.out = length(idx))
    folds[idx] <- folds_level
  }

  # optional: warn if some fold has no events
  events_per_fold <- tapply(d, folds, sum)
  if (any(events_per_fold == 0)) {
    warning("Some folds have zero events; this is unavoidable when events < nfolds.")
  }

  used <- sort(unique(folds))
  if (length(used) < nfolds) {
    warning(sprintf(
      "Requested %d folds, but only %d contained samples; reducing nfolds accordingly.",
      nfolds, length(used)
    ))
    fold_map <- setNames(seq_along(used), used)
    folds <- unname(fold_map[as.character(folds)])
  }

  attr(folds, "nfolds") <- length(unique(folds))
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

.desurv_merge_args <- function(...) {
  pieces <- list(...)
  args <- do.call(c, pieces)
  nms <- names(args)
  if (!is.null(nms)) {
    keep <- rep(TRUE, length(args))
    dup <- duplicated(nms, fromLast = TRUE)
    dup <- dup & nzchar(nms)
    keep[dup] <- FALSE
    args <- args[keep]
  }
  args
}

#' @keywords internal
.desurv_apply_validation_transform <- function(X, method, target) {
  if (method == "rank") {
    X_rank <- apply(X, 2, rank, ties.method = "average")
    dim(X_rank) <- dim(X)
    dimnames(X_rank) <- dimnames(X)
    return(X_rank)
  }
  if (method == "quant") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package `preprocessCore` is required for quantile normalization.")
    }
    target_info <- .desurv_unpack_quant_target(target)
    if (is.null(target_info$values)) {
      stop("Quantile normalization target is missing for validation data.")
    }
    gene_order <- rownames(X)
    X_norm <- preprocessCore::normalize.quantiles.use.target(
      X,
      target = target_info$values
    )
    final_genes <- if (!is.null(target_info$genes)) {
      if (length(target_info$genes) != nrow(X_norm)) {
        stop("Quantile normalization gene dimension mismatch for validation data.")
      }
      target_info$genes
    } else {
      gene_order
    }
    if (!is.null(final_genes)) {
      rownames(X_norm) <- final_genes
    }
    return(X_norm)
  }
  if (method != "none") {
    stop("Unsupported transformation specified.")
  }
  X
}

#' @keywords internal
.desurv_prepare_fold_payload <- function(
    X_full,
    y_full,
    d_full,
    dataset_full,
    idx_tr,
    idx_val,
    preprocess,
    ngene,
    genes,
    method_trans_train,
    verbose
) {
  X_tr <- X_full[, idx_tr, drop = FALSE]
  y_tr <- y_full[idx_tr]
  d_tr <- d_full[idx_tr]
  dataset_tr <- dataset_full[idx_tr]

  X_val <- X_full[, idx_val, drop = FALSE]
  y_val <- y_full[idx_val]
  d_val <- d_full[idx_val]

  if (!preprocess) {
    list(
      X_tr     = X_tr,
      y_tr     = y_tr,
      d_tr     = d_tr,
      dataset_tr = dataset_tr,
      X_val    = X_val,
      y_val    = y_val,
      d_val    = d_val,
      p_tr     = nrow(X_tr),
      n_tr     = ncol(X_tr),
      max_x_tr = {
        m <- max(X_tr)
        if (!is.finite(m) || m <= 0) 1 else m
      }
    )
  } else {
    proc <- preprocess_data(
      X                  = X_tr,
      y                  = y_tr,
      d                  = d_tr,
      dataset            = dataset_tr,
      ngene              = ngene,
      genes              = genes,
      method_trans_train = method_trans_train,
      verbose            = verbose
    )

    X_tr_proc <- proc$ex
    samp_tr   <- proc$sampInfo

    genes_keep <- proc$featInfo
    X_val_proc <- X_val[genes_keep, , drop = FALSE]
    X_val_proc[is.na(X_val_proc)] <- 0
    X_val_proc <- .desurv_apply_validation_transform(
      X_val_proc,
      method = method_trans_train,
      target = proc$transform_target
    )

    list(
      X_tr     = X_tr_proc,
      y_tr     = samp_tr$time,
      d_tr     = samp_tr$event,
      dataset_tr = samp_tr$dataset,
      X_val    = X_val_proc,
      y_val    = y_val,
      d_val    = d_val,
      p_tr     = nrow(X_tr_proc),
      n_tr     = ncol(X_tr_proc),
      max_x_tr = {
        m <- max(X_tr_proc)
        if (!is.finite(m) || m <= 0) 1 else m
      }
    )
  }
}
