#' @export
run_cv = function(X, y, delta, k, nfold, alpha, lambda = 0, eta = 0, 
                  lambdaW = 0, lambdaH = 0, fold_info, seed = 123,
                  ninit = 100, imaxit=5000, maxit = 5000, tol = 1e-6, 
                  parallel = TRUE, ncore = NULL, replace = FALSE, 
                  save = TRUE, verbose=TRUE, prefix){
  
  X = as.matrix(X)
  
  if(ncore==1){
    parallel=FALSE
  }
  
  if(parallel & is.null(ncore)){
    ncore = detectCores() - 1
  }
  
  params = set_param_grid(k=k, alpha=alpha, lambda=lambda, eta=eta, 
                          lambdaW, lambdaH, ninit=ninit, 
                          type="cv", nfold=nfold, prefix = prefix,
                          ngene = ngene, maxit=maxit, tol=tol, imaxit=imaxit)
  
  Xtrain = fold_info$Xtrain
  Xtest = fold_info$Xtest
  folds = fold_info$folds
  
  if(parallel){
    cl = parallel::makeCluster(ncore,outfile="")
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  }
  
  metrics = 
    foreach(pa=1:nrow(params), 
            .inorder = FALSE, 
            .errorhandling = 'remove', 
            .combine = 'rbind',
            .packages = c("coxNMF","survival","cvwrapr")) %dopar% 
    {#begin foreach
    
    suppressWarnings(rm(fit_cox))
    
    
    a = params$alpha[pa]
    l = params$lambda[pa]
    e = params$eta[pa]
    f = params$fold[pa]
    k = params$k[pa]
    lW = params$lambdaW[pa]
    lH = params$lambdaH[pa]
    
    Xtr = Xtrain[[f]]
    y_curr = y[folds!=f]
    delta_curr = delta[folds!=f]
    
    if(replace | !file.exists(params$file[pa])){
      print(sprintf('alpha=%f, lambda=%f, eta=%f, k=%d, f=%d\n',a,l,e,k,f))
      print("running...")
      fit_cox = run_coxNMF(X=Xtr, y=y_curr, delta=delta_curr, k=k,
                           alpha=a, lambda=l, eta=e,
                           lambdaW = lW, lambdaH = lH, seed = seed,
                           tol=tol, maxit=maxit, verbose=verbose,
                           ninit=ninit, imaxit=imaxit)
      if(save){
        save(fit_cox,file=params$file[pa])
      }
    }else{
      load(params$file[pa])
    }

    if(fit_cox$`NaN flag`){
      warning("alpha too large")
    }
    
    #compute test set metrics: loss, sloss, recon error, cindex, bic
    ytest=y[folds==f]
    dtest=delta[folds==f]
    W = fit_cox$W
    H = fit_cox$H
    beta = fit_cox$beta
    M = matrix(1,nrow=nrow(Xtest[[f]]),ncol=ncol(Xtest[[f]]))
    
    if(any(is.nan(t(Xtest[[f]]) %*% W %*% beta))){
      warning("alpha too large 2")
    }
    
    c = cvwrapr::getCindex(t(Xtest[[f]]) %*% W %*% beta, Surv(ytest, dtest))
    sl = calc_surv_loss(Xtest[[f]], M, ytest, dtest, W, beta)
    Xall=cbind(Xtrain[[f]],Xtest[[f]])
    M = matrix(1,nrow=nrow(Xall),ncol=ncol(Xall))
    slall = calc_surv_loss(Xall, M, y, delta, W, beta)
    cvl = slall - fit_cox$loss$surv_loss
    bic = -2*cvl + k*log(ncol(Xtest[[f]]))
    
    converged=fit_cox$iter < maxit
    
    data.frame(k=k,alpha=a,lambda=l,eta=e,lambdaW=lW, lambdaH=lH,
               fold=f,sloss=sl,sall=slall,iter=fit_cox$iter,strain=fit_cox$loss$surv_loss,bic=bic,c=c, cvl= cvl, converged=converged)
    
    
  }#end foreach
  
  if(parallel){
    stopCluster(cl)
  }

  return(metrics)

}
