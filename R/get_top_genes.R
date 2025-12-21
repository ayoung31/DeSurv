#' Identify top genes per latent factor
#'
#' @description
#' Compute, for each column of a factor loading matrix \code{W}, the genes
#' whose weights stand out relative to the remaining factors. Columns are
#' first scaled to unit maximum, then genes are ranked by the difference
#' between the focal column and the maximum weight achieved by any other
#' column. The top \code{ntop} genes per factor are returned alongside their
#' ranking indices and per-gene contrast scores.
#'
#' @param W Numeric matrix with genes in rows and factors in columns.
#' @param ntop Positive integer specifying how many top genes to retain per
#'   factor. Values greater than the number of genes are clamped internally.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{top_genes}: data frame of dimension \code{ntop x k} (after
#'         removing all-zero columns) containing gene identifiers per factor.
#'   \item \code{top_indices}: integer data frame mirroring \code{top_genes}
#'         but storing row indices (useful when gene names are absent).
#'   \item \code{diffs}: data frame containing the contrast scores
#'         (\code{current factor - best competing factor}) for every gene and
#'         factor considered.
#' }
#'
#' @examples
#' set.seed(1)
#' W <- matrix(rexp(50), nrow = 10)
#' rownames(W) <- paste0("gene", seq_len(10))
#' desurv_get_top_genes(W, ntop = 3)
#'
#' @export
desurv_get_top_genes <- function(W, ntop) {
  if (is.null(W) || !is.matrix(W)) {
    stop("`W` must be a non-null numeric matrix.", call. = FALSE)
  }
  if (!is.numeric(W)) {
    stop("`W` must be numeric.", call. = FALSE)
  }

  ntop <- as.integer(ntop)
  if (length(ntop) != 1L || is.na(ntop) || ntop <= 0L) {
    stop("`ntop` must be a positive integer.", call. = FALSE)
  }
  if (nrow(W) == 0L || ncol(W) == 0L) {
    stop("`W` must have at least one row and one column.", call. = FALSE)
  }
  if (ntop > nrow(W)) {
    ntop <- nrow(W)
  }

  max_per_col <- apply(W, 2, max)
  keep_cols <- max_per_col > 0
  if (!any(keep_cols)) {
    stop("All columns of `W` are zero; cannot compute top genes.", call. = FALSE)
  }
  W_use <- W[, keep_cols, drop = FALSE]
  col_names <- colnames(W_use)
  if (is.null(col_names)) {
    col_names <- paste0("factor", seq_len(ncol(W_use)))
  }

  # scale each column by its max to match the reference implementation
  scaler <- apply(W_use, 2, max)
  scaler[scaler == 0] <- 1
  W_scaled <- sweep(W_use, 2, scaler, FUN = "/")

  gene_names <- rownames(W_scaled)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(W_scaled)))
  }

  top_gene_mat <- matrix(NA_character_, nrow = ntop, ncol = ncol(W_scaled))
  top_idx_mat <- matrix(NA_integer_, nrow = ntop, ncol = ncol(W_scaled))
  diff_mat <- matrix(NA_real_, nrow = nrow(W_scaled), ncol = ncol(W_scaled))

  for (j in seq_len(ncol(W_scaled))) {
    current_col <- W_scaled[, j]
    if (sum(current_col) <= 0) {
      next
    }

    if (ncol(W_scaled) == 1L) {
      max_other <- rep(0, nrow(W_scaled))
    } else {
      other_idx <- setdiff(seq_len(ncol(W_scaled)), j)
      max_other <- apply(W_scaled[, other_idx, drop = FALSE], 1, max)
    }

    diff_vector <- current_col - max_other
    ord <- order(diff_vector, decreasing = TRUE)[seq_len(ntop)]

    top_gene_mat[, j] <- gene_names[ord]
    top_idx_mat[, j] <- ord
    diff_mat[, j] <- diff_vector
  }

  top_genes_df <- as.data.frame(top_gene_mat, stringsAsFactors = FALSE)
  names(top_genes_df) <- col_names
  rownames(top_genes_df) <- paste0("rank", seq_len(ntop))

  top_idx_df <- as.data.frame(top_idx_mat, stringsAsFactors = FALSE)
  names(top_idx_df) <- col_names
  rownames(top_idx_df) <- paste0("rank", seq_len(ntop))

  diffs_df <- as.data.frame(diff_mat, stringsAsFactors = FALSE)
  names(diffs_df) <- col_names
  rownames(diffs_df) <- gene_names

  list(
    top_genes = top_genes_df,
    top_indices = top_idx_df,
    diffs = diffs_df
  )
}

.desurv_normalize_ntop <- function(ntop, max_genes = NULL) {
  if (is.null(ntop) || length(ntop) == 0L) {
    return(NULL)
  }
  ntop_num <- suppressWarnings(as.numeric(ntop)[1L])
  if (is.na(ntop_num) || ntop_num <= 0) {
    return(NULL)
  }
  ntop_int <- as.integer(ntop_num)
  if (!is.null(max_genes)) {
    ntop_int <- min(ntop_int, max_genes)
  }
  if (ntop_int <= 0L) {
    NULL
  } else {
    ntop_int
  }
}

.desurv_lp_with_top_genes <- function(fit, X_new, ntop) {
  ntop_use <- .desurv_normalize_ntop(ntop, nrow(fit$W))
  if (is.null(ntop_use)) {
    return(predict(fit, newdata = X_new, type = "lp"))
  }

  top_info <- desurv_get_top_genes(fit$W, ntop_use)
  idx_mat <- top_info$top_indices
  if (is.null(idx_mat) || !nrow(idx_mat)) {
    return(predict(fit, newdata = X_new, type = "lp"))
  }

  idx <- unique(as.integer(unlist(idx_mat, use.names = FALSE)))
  idx <- idx[!is.na(idx) & idx >= 1L & idx <= nrow(fit$W)]
  if (!length(idx)) {
    return(predict(fit, newdata = X_new, type = "lp"))
  }

  theta = fit$W %*% fit$beta

  theta_sub = theta[idx, ,drop=FALSE]
  theta_norm = sqrt(sum(theta_sub^2))
  if (!is.finite(theta_norm)) {
    theta_sub[] = 0
  } else if (theta_norm > 0) {
    theta_sub = theta_sub / theta_norm
  }

  X_sub <- X_new[idx, , drop = FALSE]
  drop(t(X_sub) %*% theta_sub)
}
