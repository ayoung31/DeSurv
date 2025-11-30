#' Predictions for fitted DeSurv models
#'
#' @param object A `desurv_fit` object returned by [desurv_fit()].
#' @param newdata Optional numeric matrix, data frame, or `desurv_data`
#'   providing expression profiles for which predictions should be
#'   generated. When `NULL`, the training matrix stored in `object$data`
#'   is used.
#' @param type Character string indicating the desired prediction scale.
#'   Supported values are `"lp"` (linear predictor, default) and `"risk"`
#'   (exponentiated linear predictor).
#' @param ... Unused; included for S3 completeness.
#'
#' @return A numeric vector with one element per column of `newdata`.
#' @export
predict.desurv_fit <- function(object,
                               newdata = NULL,
                               type = c("lp", "risk"),
                               ...) {
  stopifnot(inherits(object, "desurv_fit"))
  type <- match.arg(type)

  X <- .desurv_prepare_newdata(newdata, object)
  Z <- t(X) %*% object$W
  lp <- drop(Z %*% object$beta)

  switch(
    type,
    lp   = lp,
    risk = exp(lp)
  )
}

#' @rdname predict.desurv_fit
#' @export
coef.desurv_fit <- function(object, ...) {
  stopifnot(inherits(object, "desurv_fit"))
  object$beta
}

.desurv_prepare_newdata <- function(newdata, object) {
  template <- object$data$X
  preprocess_meta <- object$preprocess

  if (is.null(newdata)) {
    return(template)
  }

  if (inherits(newdata, "desurv_data")) {
    return(newdata$X)
  }

  X <- as.matrix(newdata)
  if (!is.numeric(X)) {
    stop("`newdata` must be coercible to a numeric matrix.")
  }

  train_genes <- if (!is.null(preprocess_meta$genes)) {
    preprocess_meta$genes
  } else {
    rownames(template)
  }

  if (!is.null(train_genes)) {
    if (is.null(rownames(X))) {
      if (nrow(X) != length(train_genes)) {
        stop("`newdata` must have rownames matching the training features.", call. = FALSE)
      }
      rownames(X) <- train_genes
    }
    missing_genes <- setdiff(train_genes, rownames(X))
    if (length(missing_genes)) {
      stop(
        "The following required genes are missing from `newdata`: ",
        paste(head(missing_genes, 10L), collapse = ", "),
        if (length(missing_genes) > 10) ", ..."
      )
    }
    X <- X[train_genes, , drop = FALSE]
  } else if (nrow(X) != nrow(template)) {
    stop("`newdata` must have ", nrow(template), " rows to match the training data.")
  }

  if (!is.null(preprocess_meta)) {
    X[is.na(X)] <- 0
    method <- preprocess_meta$method_trans_train
    if (isTRUE(method == "rank")) {
      X_rank <- apply(X, 2, rank, ties.method = "average")
      dim(X_rank) <- dim(X)
      dimnames(X_rank) <- dimnames(X)
      X <- X_rank
    } else if (isTRUE(method == "quant")) {
      if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        stop("Package `preprocessCore` is required for quantile normalization.")
      }
      target_info <- .desurv_unpack_quant_target(
        preprocess_meta$transform_target,
        fallback_genes = train_genes
      )
      if (is.null(target_info$values)) {
        stop("Quantile normalization target is missing from the fitted model.", call. = FALSE)
      }
      gene_order <- rownames(X)
      X <- preprocessCore::normalize.quantiles.use.target(
        X,
        target = target_info$values
      )
      final_genes <- if (!is.null(target_info$genes)) {
        if (nrow(X) != length(target_info$genes)) {
          stop("Quantile normalization gene dimension mismatch at prediction time.", call. = FALSE)
        }
        target_info$genes
      } else {
        gene_order
      }
      if (!is.null(final_genes)) {
        rownames(X) <- final_genes
      }
    }
  }

  storage.mode(X) <- "double"
  X
}
