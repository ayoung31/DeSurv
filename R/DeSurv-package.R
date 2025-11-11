#' DeSurv: Survival-Driven Nonnegative Matrix Factorization
#'
#' @description
#' DeSurv couples nonnegative matrix factorisation with a Cox proportional
#' hazards model so that survival information steers the latent factors that
#' summarise high-dimensional molecular data. The package exposes
#' multi-start optimisers, alpha warm-start paths, hyperparameter
#' cross-validation helpers, prediction utilities, and validation helpers for
#' expression matrices typically encountered in RNA-seq studies.
#'
#' @docType package
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib DeSurv, .registration=TRUE
#' @importFrom Rcpp evalCpp
## usethis namespace: end
NULL
