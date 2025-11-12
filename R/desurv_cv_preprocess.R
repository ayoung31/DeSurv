#' Preprocess expression data then run cross-validated DeSurv (deprecated)
#'
#' @description
#' This convenience wrapper has been superseded by the
#'   \code{preprocess = TRUE} option inside \code{\link{desurv_cv}}. It
#'   remains for backward compatibility but will emit a soft deprecation
#'   warning.
#'
#' @inheritParams preprocess_data
#' @inheritParams desurv_cv
#' @param ... Additional arguments passed to \code{\link{desurv_cv}} (e.g.,
#'   \code{k_grid}, \code{alpha_grid}, \code{lambda_grid}, etc.).
#'
#' @return A \code{"desurv_cv"} object identical to calling
#'   \code{desurv_cv(preprocess = TRUE, ...)}.
#'
#' @seealso \code{\link{desurv_cv}}, \code{\link{preprocess_data}}
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
  .Deprecated("desurv_cv(preprocess = TRUE, ...)")

  desurv_cv(
    X                  = X,
    y                  = y,
    d                  = d,
    dataset            = dataset,
    samp_keeps         = samp_keeps,
    ngene              = ngene,
    genes              = genes,
    method_trans_train = method_trans_train,
    preprocess         = TRUE,
    ...
  )
}
