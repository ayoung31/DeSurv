#' Select highly expressed and variable genes
#'
#' @param X Numeric matrix with genes as rows and samples as columns.
#' @param ngene Positive integer indicating how many genes to retain.
#'
#' @return A numeric matrix containing the top `ngene` genes ranked by
#'   average expression and variability (ties broken by combined rank).
#' @export
gene_filter <- function(X, ngene = 1000) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be a numeric matrix.")
  if (is.null(rownames(X))) stop("`X` must have gene identifiers as rownames.")

  ngene <- as.integer(ngene)
  if (length(ngene) != 1L || is.na(ngene) || ngene <= 0L) {
    stop("`ngene` must be a positive integer.")
  }

  ngene <- min(ngene, nrow(X))

  mean_expr <- rowMeans(X, na.rm = TRUE)
  var_expr  <- apply(X, 1L, stats::var, na.rm = TRUE)

  expr_rank <- rank(-mean_expr, ties.method = "average")
  var_rank  <- rank(-var_expr,  ties.method = "average")
  combined  <- expr_rank + var_rank

  keep_idx <- order(combined)[seq_len(ngene)]
  X[keep_idx, , drop = FALSE]
}
