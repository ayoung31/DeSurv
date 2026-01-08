.desurv_make_quant_target <- function(values, genes) {
  if (is.null(values)) {
    return(NULL)
  }
  values_num <- as.numeric(values)
  list(
    values = values_num,
    genes  = if (is.null(genes)) NULL else as.character(genes)
  )
}

.desurv_unpack_quant_target <- function(target, fallback_genes = NULL) {
  if (is.null(target)) {
    return(list(values = NULL, genes = fallback_genes))
  }
  if (is.list(target)) {
    values <- target$values
    if (is.null(values) && !is.null(target$target)) {
      values <- target$target
    }
    genes <- target$genes
    if (is.null(genes)) {
      genes <- fallback_genes
    }
    return(list(
      values = if (is.null(values)) NULL else as.numeric(values),
      genes  = if (is.null(genes)) NULL else as.character(genes)
    ))
  }
  list(
    values = as.numeric(target),
    genes  = fallback_genes
  )
}

#' Preprocess multi-cohort expression data for DeSurv models
#'
#' @description
#' Filters samples, restricts genes (either via `gene_filter()` per dataset or
#' via a supplied gene list), and applies an optional within-sample
#' transformation prior to calling downstream DeSurv routines.
#'
#' @param X Numeric matrix (genes x samples) of expression values.
#' @param y Numeric vector of survival times, length `ncol(X)`.
#' @param d Numeric (0/1) event indicator vector, length `ncol(X)`.
#' @param dataset Vector identifying the cohort/dataset for each sample.
#' @param samp_keeps Optional logical/integer/character index of samples to
#'   retain before preprocessing. Defaults to all samples.
#' @param ngene Integer number of genes to retain per dataset when `genes` is
#'   `NULL`.
#' @param genes Optional character vector of genes to keep (overrides `ngene`
#'   filtering).
#' @param method_trans_train Character string: `"rank"`, `"quant"`, or `"none"`
#'   describing the per-sample transformation applied after filtering.
#' @param verbose Logical; if `TRUE`, emit preprocessing progress messages.
#'
#' @return A list with components
#'   \describe{
#'     \item{ex}{Processed expression matrix.}
#'     \item{sampInfo}{Data frame containing `time`, `event`, and `dataset`.}
#'     \item{featInfo}{Character vector of retained gene identifiers.}
#'     \item{samp_keeps}{Integer indices of samples retained.}
#'   }
#'
#' @export
preprocess_data <- function(
    X,
    y,
    d,
    dataset,
    samp_keeps         = NULL,
    ngene              = 1000,
    genes              = NULL,
    method_trans_train = c("rank", "quant", "none"),
    verbose            = TRUE
) {
  method_trans_train <- match.arg(method_trans_train)
  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("`X` must be a numeric matrix.")
  }

  n_samples <- ncol(X)
  if (length(y) != n_samples || length(d) != n_samples || length(dataset) != n_samples) {
    stop("`X`, `y`, `d`, and `dataset` must have matching sample counts.")
  }


  # Finding 1 fix: Validate y, d, and dataset
  # Validate y: must be numeric, finite, positive, no NA
  if (!is.numeric(y)) {
    stop("`y` must be a numeric vector of survival times.", call. = FALSE)
  }
  if (anyNA(y)) {
    stop("`y` cannot contain NA values.", call. = FALSE)
  }
  if (!all(is.finite(y))) {
    stop("`y` cannot contain Inf or NaN values.", call. = FALSE)
  }
  if (any(y <= 0)) {
    stop("`y` (survival times) must be positive (> 0).", call. = FALSE)
  }

  # Validate d: must be numeric or integer, values in {0, 1}, no NA
  if (!is.numeric(d) && !is.integer(d)) {
    stop("`d` must be a numeric or integer event indicator vector.", call. = FALSE)
  }
  if (anyNA(d)) {
    stop("`d` cannot contain NA values.", call. = FALSE)
  }
  d_unique <- unique(d)
  if (!all(d_unique %in% c(0, 1))) {
    stop("`d` must contain only 0 (censored) and 1 (event) values.", call. = FALSE)
  }

  # Validate dataset: no NA
  if (anyNA(dataset)) {
    stop("`dataset` cannot contain NA values.", call. = FALSE)
  }

  # Normalize sample subset indices
  normalize_samp_idx <- function(idx, n, sample_names) {
    if (is.logical(idx)) {
      if (length(idx) != n) stop("Logical `samp_keeps` must match number of samples.")
      which(idx)
    } else if (is.numeric(idx)) {
      idx <- as.integer(idx)
      if (anyNA(idx) || any(idx < 1L | idx > n)) stop("Numeric `samp_keeps` out of bounds.")
      unique(idx)
    } else if (is.character(idx)) {
      if (is.null(sample_names)) stop("Character `samp_keeps` requires column names on `X`.")
      pos <- match(idx, sample_names)
      if (anyNA(pos)) stop("Character `samp_keeps` contains unknown sample IDs.")
      unique(pos)
    } else {
      stop("`samp_keeps` must be logical, numeric, character, or NULL.")
    }
  }

  keep_idx <- if (is.null(samp_keeps)) {
    seq_len(n_samples)
  } else {
    normalize_samp_idx(samp_keeps, n_samples, colnames(X))
  }

  if (length(keep_idx) == 0L) {
    stop("No samples selected after applying `samp_keeps`.")
  }

  X <- X[, keep_idx, drop = FALSE]
  y <- y[keep_idx]
  d <- d[keep_idx]
  dataset <- dataset[keep_idx]
  sample_ids <- colnames(X)

  # Restrict genes
  if (is.null(genes)) {
    datasets <- unique(dataset)
    per_dataset <- lapply(
      datasets,
      function(ds) {
        cols <- dataset == ds
        gene_filter(X = X[, cols, drop = FALSE], ngene = ngene)
      }
    )
    keep_genes <- Reduce(intersect, lapply(per_dataset, rownames))
    if (length(keep_genes) == 0L) {
      stop("No overlapping genes found across datasets after filtering.")
    }
    if (isTRUE(verbose)) {
      message(
        paste0(
          "After preprocessing and dataset merging, ",
          length(keep_genes),
          " genes were kept"
        )
      )
    }
    X <- X[keep_genes, , drop = FALSE]
  } else {
    genes <- unique(genes)
    keep_genes <- genes[genes %in% rownames(X)]
    if (length(keep_genes) == 0L) {
      stop("None of the supplied `genes` are present in `X`.")
    }
    # Finding 6 fix: Warn when explicit genes are dropped
    dropped_genes <- setdiff(genes, keep_genes)
    if (length(dropped_genes) > 0L) {
      warning(
        sprintf(
          "%d of %d requested genes not found in X and were dropped: %s%s",
          length(dropped_genes),
          length(genes),
          paste(head(dropped_genes, 5L), collapse = ", "),
          if (length(dropped_genes) > 5L) ", ..." else ""
        ),
        call. = FALSE
      )
    }
    X <- X[keep_genes, , drop = FALSE]
  }

  # Finding 2 fix: Replace all non-finite values (NA, NaN, Inf, -Inf) with 0
  # This prevents corrupted variance calculations and quantile targets
  non_finite_mask <- !is.finite(X)
  if (any(non_finite_mask)) {
    n_nonfinite <- sum(non_finite_mask)
    if (isTRUE(verbose)) {
      message(sprintf("Replaced %d non-finite values (NA/NaN/Inf) with 0.", n_nonfinite))
    }
    X[non_finite_mask] <- 0
  }

  transform_target <- NULL
  if (method_trans_train == "rank") {
    X_rank <- apply(X, 2, rank, ties.method = "average")
    dim(X_rank) <- dim(X)
    dimnames(X_rank) <- dimnames(X)
    X <- X_rank
  } else if (method_trans_train == "quant") {
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package `preprocessCore` is required for quantile normalization.")
    }
    gene_order <- rownames(X)
    # Store the training target and gene order for prediction-time reuse.
    target_values <- preprocessCore::normalize.quantiles.determine.target(X)
    transform_target <- .desurv_make_quant_target(target_values, gene_order)
    X <- preprocessCore::normalize.quantiles.use.target(
      X,
      target = transform_target$values
    )
    if (!is.null(transform_target$genes)) {
      if (nrow(X) != length(transform_target$genes)) {
        stop("Quantile normalization changed gene dimensions unexpectedly.")
      }
      rownames(X) <- transform_target$genes
    }
  } else if (method_trans_train != "none") {
    stop("Unsupported transformation specified.")
  }

  sample_ids <- colnames(X)

  data.frame_samp <- data.frame(
    time    = y,
    event   = d,
    dataset = dataset,
    row.names = sample_ids,
    check.names = FALSE
  )

  list(
    ex        = X,
    sampInfo          = data.frame_samp,
    featInfo          = rownames(X),
    samp_keeps        = keep_idx,
    transform_target  = transform_target,
    method_trans_train = method_trans_train
  )
}
