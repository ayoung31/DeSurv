cv.Desurv <- function(
    X, y, delta,
    params, ninit, alpha, tol, maxit,
    verbose = TRUE,
    path_fits
){
  k                  = params$k
  lambda             = params$lambda
  eta                = params$eta
  lambdaW            = params$lambdaW
  lambdaH            = params$lambdaH
  seeds              = 1:ninit
  flag_exist         = FALSE
  f                  = params$fold

  dir.create(dirname(path_fits), recursive = TRUE, showWarnings = FALSE)
  # if(file.exists(path_fits)){
  #   bundle_old = readRDS(path_fits)
  #   exists = unlist(lapply(bundle_old,function(x) names(x$fits)==as.character(alpha)))
  #   if(all(exists)){
  #     if(length(bundle_old)==ninit){
  #       return(path_fits)
  #     }else{
  #       seeds = setdiff(1:ninit,1:length(bundle_old))
  #       flag_exist = TRUE
  #     }
  #   }
  # 
  # }
  
  # cl <- parallel::makeCluster(NINIT,setup_strategy="sequential")
  # registerDoParallel(cl)
  doMC::registerDoMC(cores = ninit) 

  bundle <- foreach::foreach(
    i = seeds,
    .combine  = 'c',
    .inorder  = FALSE,
    .export   = c("run_coxNMF", "alpha", "tol", "maxit"),
    .errorhandling = "remove",
    .packages = c("coxNMF","survival")  # add packages run_coxNMF needs, e.g. "Matrix", "survival"
  ) %dopar% {
    fit_prev = NULL
    fits = list()
    for (a in alpha) {
      print(sprintf("seed %d, alpha %.2f",i,a))
      if (is.null(fit_prev)) {
        fit <- run_coxNMF(
          X=X, y=y, delta=delta, k=k,
          alpha=a, lambda=lambda, eta=eta,
          lambdaW=lambdaW, lambdaH=lambdaH,
          seed=i, tol=tol, maxit=maxit, verbose=verbose,
          ninit=1, imaxit=maxit
        )
      } else {
        fit <- run_coxNMF(
          X=X, y=y, delta=delta, k=k,
          alpha=a, lambda=lambda, eta=eta,
          lambdaW=lambdaW, lambdaH=lambdaH,
          tol=tol, maxit=maxit, verbose=verbose,
          ninit=1, imaxit=maxit,
          W0=fit_prev$W, H0=fit_prev$H, beta0=fit_prev$beta
        )
      }
      
      fits[[as.character(a)]] <- fit
      fit_prev <- fit
    }
  
    
    # assemble / merge bundle (unchanged, plus best_seed in meta)
    meta <- list(
      k=k, lambda=lambda, eta=eta, lambdaW=lambdaW, lambdaH=lambdaH,
      alpha_vec = alpha, maxit=maxit, tol=tol,
      ngene = nrow(X), nsamples = ncol(X), fold=f,
      seed = i
    )
  
    setNames(list(list(meta = meta, fits = fits)), as.character(i))
  }
  
  # stopCluster(cl)
  
  # if(flag_exist){
  #   bundle = c(bundle_old,bundle)
  # }
  
  
  tmp <- paste0(path_fits, ".tmp")
  saveRDS(bundle, tmp, compress = "xz")
  file.rename(tmp, path_fits)
  
  path_fits
}
