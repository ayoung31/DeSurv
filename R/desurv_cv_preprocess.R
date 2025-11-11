#' Preprocess expression data and run cross-validated DeSurv
#'
#' @description
#' Convenience wrapper that first calls \code{\link{preprocess_data}} on raw
#' multi-cohort expression inputs and then passes the processed expression
#' matrix and survival outcomes to \code{\link{desurv_cv}}.
#'
#' @param X Numeric matrix (genes x samples) of expression values.
#' @param y Numeric vector of survival/event times aligned with the columns of
#'   \code{X}.
#' @param d Numeric (0/1) event indicator vector aligned with \code{X}.
#' @param dataset Vector identifying the cohort/dataset for each sample (length
#'   must equal \code{ncol(X)}).
#' @param samp_keeps Optional logical/integer/character index of samples to
#'   retain prior to preprocessing.
#' @param ngene Integer; if \code{genes} is \code{NULL}, the number of
#'   highly expressed/variable genes to keep per dataset (see
#'   \code{\link{preprocess_data}}).
#' @param genes Optional character vector of gene IDs to restrict to. If
#'   \code{NULL}, top \code{ngene} genes are selected by
#'   \code{\link{gene_filter}} within each dataset.
#' @param method_trans_train Character; one of \code{"rank"},
#'   \code{"quant"}, or \code{"none"}, passed to \code{preprocess_data()}
#'   to control the within-sample transformation.
#' @param ... Additional arguments passed directly to \code{\link{desurv_cv}}
#'   (e.g., \code{k_grid}, \code{alpha_grid}, \code{lambda_grid},
#'   \code{nu_grid}, \code{lambdaW_grid}, \code{lambdaH_grid},
#'   \code{n_starts}, \code{nfolds}, \code{rule}, etc.).
#'
#' @return
#' An object of class \code{"desurv_cv"} as returned by \code{desurv_cv()},
#' with an additional component \code{$preprocess} containing
#' preprocessing metadata (selected genes, transformation, and sample
#' subset).
#'
#' @seealso
#' \code{\link{preprocess_data}}, \code{\link{desurv_cv}}
#'
#' @export
desurv_cv_preprocess <- function(
    X,
    y,
    d,
    dataset,
    samp_keeps         = NULL,
    ngene              = 1000,
    genes              = NULL,
    method_trans_train = c("rank", "quant", "none"),
    ...
) {
  method_trans_train <- match.arg(method_trans_train)

  # Run preprocessing (sample/gene filtering + transformation)
  data_proc <- preprocess_data(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    samp_keeps         = samp_keeps,
    ngene              = ngene,
    genes              = genes,
    method_trans_train = method_trans_train
  )

  # Check survival columns
  required_cols <- c("time", "event")
  if (!all(required_cols %in% colnames(data_proc$sampInfo))) {
    stop(
      "Expected columns 'time' and 'event' in data$sampInfo after preprocessing."
    )
  }

  # Extract processed X and aligned y/d from sampInfo
  X <- data_proc$ex
  y <- data_proc$sampInfo$time
  d <- data_proc$sampInfo$event

  # Call the main CV engine on the processed data
  cv_fit <- desurv_cv(
    X = X,
    y = y,
    d = d,
    ...
  )

  # Attach preprocessing metadata for later inspection
  cv_fit$preprocess <- list(
    ngene              = ngene,
    genes              = rownames(X),
    method_trans_train = method_trans_train,
    samp_keeps         = data_proc$samp_keeps
  )

  cv_fit
}
