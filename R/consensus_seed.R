#' Build a consensus initialization from multiple DeSurv fits
#'
#' @description
#' Aggregate the top genes that appear in each latent factor across several
#' independently initialized \code{desurv_fit} objects, compute a consensus
#' co-occurrence matrix, and derive a deterministic \eqn{W_0} (and matching
#' \eqn{H_0}) that can be supplied to \code{\link{desurv_fit}} for a final,
#' fully optimized run. The consensus weights are obtained by clustering the
#' co-occurrence matrix, normalising the resulting cluster profiles, and
#' solving non-negative least squares problems to align \eqn{H_0} with the
#' observed data matrix \code{X}.
#'
#' @param fits A list of \code{desurv_fit} objects (or a single fit) obtained
#'   from running \code{desurv_fit()} with different random seeds.
#' @param X Numeric matrix containing the training data used for each fit
#'   (genes in rows, samples in columns). Row order should match the
#'   \code{W} matrices stored in \code{fits}. If rownames are available they
#'   will be used to align \code{X} to the consensus initialization.
#' @param ntop Positive integer specifying how many top genes per factor are
#'   considered when computing the co-occurrence matrix.
#' @param k Optional integer overriding the rank taken from the supplied fits.
#'   Defaults to the number of columns in the first fit's \code{W}.
#' @param clustering Character string passed to \code{\link[stats]{hclust}}
#'   to control the agglomerative linkage (default \code{"average"}).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{W0}: consensus basis matrix (genes x factors).
#'   \item \code{H0}: non-negative coefficients computed via NNLS so that
#'         \code{W0 \%*\% H0} approximates \code{X}.
#'   \item \code{beta0}: zero vector of length \code{k} (for convenience when
#'         seeding \code{desurv_fit()}).
#'   \item \code{consensus}: normalised gene co-occurrence matrix.
#'   \item \code{counts}: raw co-occurrence counts.
#'   \item \code{frequency}: per-gene selection counts (diagonal of
#'         \code{counts}).
#'   \item \code{clusters}: integer cluster labels used to build \code{W0}.
#'   \item \code{ntop}: number of top genes used per factor.
#'   \item \code{n_factors}: total number of short-run factors contributing to
#'         the consensus.
#' }
#'
#' @importFrom stats as.dist cutree hclust
#' @export
desurv_consensus_seed <- function(fits,
                                  X,
                                  ntop,
                                  k = NULL,
                                  clustering = c("average", "complete", "single")) {
  if (inherits(fits, "desurv_fit")) {
    fits <- list(fits)
  }
  if (!is.list(fits) || !length(fits)) {
    stop("`fits` must be a non-empty list of `desurv_fit` objects.", call. = FALSE)
  }
  if (!all(vapply(fits, inherits, logical(1), what = "desurv_fit"))) {
    stop("All elements of `fits` must inherit from `desurv_fit`.", call. = FALSE)
  }

  clustering <- match.arg(clustering)
  ntop <- as.integer(ntop)
  if (length(ntop) != 1L || is.na(ntop) || ntop <= 0L) {
    stop("`ntop` must be a positive integer.", call. = FALSE)
  }

  template_W <- fits[[1]]$W
  if (is.null(template_W) || !is.matrix(template_W)) {
    stop("Each fit must store a numeric `W` matrix.", call. = FALSE)
  }
  gene_names <- rownames(template_W)
  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(nrow(template_W)))
  }
  # ensure all fits share the same row order so that integer indices align
  same_layout <- vapply(
    fits,
    function(fit) {
      W <- fit$W
      if (is.null(rownames(W))) {
        identical(dim(W), dim(template_W))
      } else {
        identical(rownames(W), rownames(template_W))
      }
    },
    logical(1)
  )
  if (!all(same_layout)) {
    stop("All fits must have the same gene ordering (matching rownames(W)).", call. = FALSE)
  }

  k_fit <- ncol(template_W)
  if (is.null(k)) {
    k <- k_fit
  } else {
    k <- as.integer(k)
    if (length(k) != 1L || is.na(k) || k <= 0L) {
      stop("`k` must be a positive integer when supplied.", call. = FALSE)
    }
    if (k > k_fit) {
      warning("Requested k exceeds the rank in the fits; using k from fits instead.")
      k <- k_fit
    }
  }
  if (k > nrow(template_W)) {
    warning("Requested k exceeds the number of genes; clamping to nrow(W).")
    k <- nrow(template_W)
  }

  consensus_info <- .desurv_compute_consensus_matrix(fits, ntop, gene_names)
  consensus <- consensus_info$consensus
  counts <- consensus_info$counts
  clusters <- .desurv_cluster_consensus(consensus, k, method = clustering)
  W0 <- .desurv_build_W0(consensus, clusters)

  X_aligned <- .desurv_align_matrix_rows(X, rownames(W0))
  H0 <- .desurv_compute_nnls_H0(W0, X_aligned)
  beta0 <- rep(0, ncol(W0))

  list(
    W0 = W0,
    H0 = H0,
    beta0 = beta0,
    consensus = consensus,
    counts = counts,
    frequency = diag(counts),
    clusters = clusters,
    ntop = ntop,
    n_factors = consensus_info$n_factors
  )
}

.desurv_compute_consensus_matrix <- function(fits, ntop, gene_names) {
  p <- nrow(fits[[1]]$W)
  ntop <- min(ntop, p)
  counts <- matrix(0, nrow = p, ncol = p)
  rownames(counts) <- colnames(counts) <- gene_names
  total_sets <- 0L

  for (fit in fits) {
    top_idx_df <- desurv_get_top_genes(fit$W, ntop)$top_indices
    for (factor_col in seq_len(ncol(top_idx_df))) {
      idx <- unique(as.integer(top_idx_df[, factor_col]))
      idx <- idx[!is.na(idx) & idx >= 1L & idx <= p]
      if (!length(idx)) {
        next
      }
      counts[idx, idx] <- counts[idx, idx] + 1
      total_sets <- total_sets + 1L
    }
  }

  if (total_sets == 0L) {
    stop("Failed to collect any top genes from the supplied fits.", call. = FALSE)
  }

  consensus <- counts / total_sets
  list(consensus = consensus, counts = counts, n_factors = total_sets)
}

.desurv_cluster_consensus <- function(consensus, k, method = "average") {
  p <- nrow(consensus)
  if (p == 1L) {
    return(rep(1L, 1))
  }

  k <- max(1L, min(as.integer(k), p))
  dist_full <- 1 - consensus
  dist_full[dist_full < 0] <- 0
  diag(dist_full) <- 0
  dist_mat <- stats::as.dist(dist_full)
  if (all(dist_mat == 0)) {
    clusters <- rep(seq_len(k), length.out = p)
  } else {
    hc <- stats::hclust(dist_mat, method = method)
    clusters <- stats::cutree(hc, k = k)
  }
  clusters
}

.desurv_build_W0 <- function(consensus, clusters) {
  uniq_clusters <- sort(unique(clusters))
  p <- nrow(consensus)
  k <- length(uniq_clusters)
  W0 <- matrix(0, nrow = p, ncol = k,
               dimnames = list(rownames(consensus),
                               paste0("consensus_", seq_len(k))))

  for (j in seq_along(uniq_clusters)) {
    cl <- uniq_clusters[j]
    idx <- which(clusters == cl)
    profile <- if (length(idx) == 1L) {
      consensus[idx, , drop = FALSE]
    } else {
      colMeans(consensus[idx, , drop = FALSE])
    }
    W0[, j] <- as.numeric(profile)
  }

  col_max <- apply(W0, 2, max)
  col_max[col_max == 0] <- 1
  W0 <- sweep(W0, 2, col_max, "/", check.margin = FALSE)
  W0
}

.desurv_align_matrix_rows <- function(X, target_names) {
  X <- as.matrix(X)
  if (is.null(target_names)) {
    return(X)
  }
  if (is.null(rownames(X))) {
    if (nrow(X) != length(target_names)) {
      stop("Row count of `X` does not match the consensus gene set.", call. = FALSE)
    }
    rownames(X) <- target_names
    return(X)
  }
  idx <- match(target_names, rownames(X))
  if (anyNA(idx)) {
    stop("`X` is missing rows required by the consensus gene set.", call. = FALSE)
  }
  X[idx, , drop = FALSE]
}

.desurv_compute_nnls_H0 <- function(W0, X) {
  W0 <- as.matrix(W0)
  X <- as.matrix(X)
  if (!is.numeric(W0) || !is.numeric(X)) {
    stop("`W0` and `X` must be numeric matrices.", call. = FALSE)
  }
  if (nrow(W0) != nrow(X)) {
    stop("`W0` and `X` must have the same number of rows.", call. = FALSE)
  }
  if (any(W0 < 0) || any(X < 0)) {
    warning("Negative entries detected in W0 or X; NNLS may clip them to zero.")
  }
  k <- ncol(W0)
  n <- ncol(X)
  H0 <- matrix(0, nrow = k, ncol = n,
               dimnames = list(colnames(W0), colnames(X)))

  residual <- NULL
  for (j in seq_len(n)) {
    nnls_fit <- .desurv_nnls(W0, X[, j])
    residual <- nnls_fit$residual
    H0[, j] <- nnls_fit$x
  }
  H0
}

.desurv_nnls <- function(A, b, tol = 1e-10, maxiter = 5000) {
  m <- nrow(A)
  n <- ncol(A)
  if (length(b) != m) {
    stop("Length of vector `b` must match nrow(A).", call. = FALSE)
  }
  passive <- rep(FALSE, n)
  x <- numeric(n)
  residual <- b
  w <- t(A) %*% residual
  iter <- 0L

  while (any(!passive & w > tol) && iter < maxiter) {
    t_idx <- which.max(ifelse(passive, -Inf, w))
    passive[t_idx] <- TRUE
    repeat {
      Ap <- A[, passive, drop = FALSE]
      if (!ncol(Ap)) {
        break
      }
      z <- numeric(n)
      z[passive] <- tryCatch(
        qr.coef(qr(Ap), b),
        error = function(e) rep(0, ncol(Ap))
      )
      if (all(z[passive] >= -tol)) {
        x <- z
        break
      }
      negative <- passive & (z < 0)
      if (!any(negative)) {
        x <- z
        break
      }
      alpha <- min(x[negative] / (x[negative] - z[negative]))
      alpha <- if (is.finite(alpha)) alpha else 0
      x <- x + alpha * (z - x)
      to_remove <- passive & (x < tol)
      passive[to_remove] <- FALSE
      x[to_remove] <- 0
      if (!any(passive)) {
        break
      }
    }
    residual <- b - A %*% x
    w <- t(A) %*% residual
    iter <- iter + 1L
  }
  if (iter >= maxiter) {
    warning("NNLS reached maximum iterations; results may be approximate.")
  }
  list(x = pmax(x, 0), residual = residual)
}
