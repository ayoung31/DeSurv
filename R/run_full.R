#' @export
run_full = function(X, y, delta, k, alpha, lambda = 0, eta = 0,
                    lambdaW = 0, lambdaH = 0,
                    ninit = 100, imaxit = 30, maxit = 3000, tol = 1e-5, 
                    parallel = TRUE, ncore = NULL, replace = FALSE, 
                    save = TRUE, verbose=TRUE, prefix,test=FALSE){
  
  X=as.matrix(X)
  
  if(ncore==1){
    parallel=FALSE
  }
  
  if(parallel & is.null(ncore)){
    ncore = detectCores() - 1
  }
  
  params = set_param_grid(k=k, alpha=alpha, lambda=lambda, eta=eta, 
                          lambdaW, lambdaH, ninit=ninit, 
                          type="full", prefix=prefix,
                          ngene = ngene, maxit=maxit, tol=tol, imaxit=imaxit)
  
  if(parallel){
    cl = parallel::makeCluster(ncore,outfile="")
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  }
  
  if(test){
    comb = "list"
  }else{
    comb = "rbind"
  }
  
  metrics = 
    foreach(pa=1:nrow(params), .inorder = FALSE, .errorhandling = 'remove', 
            .combine = comb, .packages = c("coxNMF","survival","cvwrapr")) %dopar% {
    
    a = params$alpha[pa]
    l = params$lambda[pa]
    e = params$eta[pa]
    k = params$k[pa]
    lW = params$lambdaW[pa]
    lH = params$lambdaH[pa]
    
    
    
    if(replace | !file.exists(params$file[pa])){
      print(sprintf('alpha=%f, lambda=%f, eta=%f, k=%d\n',a,l,e,k))
      print("running...")
      fit_cox = run_coxNMF(X=X, y=y, delta=delta, k=k, 
                           alpha=a, lambda=l, eta=e, 
                           lambdaW = lW, lambdaH = lH,
                           tol=tol, maxit=maxit, verbose=verbose,
                           ninit=ninit, imaxit=imaxit)
      if(save){
        save(fit_cox,file=params$file[pa])
      }
    }else{
      #print("loading...")
      load(params$file[pa])
    }
    
    if(fit_cox$`NaN flag`){
      warning("alpha too large")
    }

    #primary metrics to output
    M=matrix(1,nrow=nrow(X),ncol=ncol(X))
    W = fit_cox$W
    H = fit_cox$H
    beta = fit_cox$beta
    
    if(any(is.nan(t(X) %*% W %*% beta))){
      warning("alpha too large 2")
    }
    
    c = cvwrapr::getCindex(t(X) %*% W %*% beta, Surv(y, delta))
    loss = fit_cox$loss
    ol = loss$loss
    sl = loss$surv_loss
    sl_std=loss$surv_loss/loss$std_surv
    nl = loss$nmf_loss
    nl_std=loss$nmf_loss/loss$std_nmf
    bic = -2*sl + k*log(ncol(X))
    converged=fit_cox$iter < maxit
    
    mets=data.frame(k=k,alpha=a,lambda=l,eta=e, lambdaW=lW, lambdaH=lH, c=c,loss=ol,
                    sloss=sl,sloss_std=sl_std,nloss=nl,nloss_std=nl_std,
                    penW=loss$penalty_W,penH=loss$penalty_H,pen_beta=loss$penalty_beta,
                    bic=bic,converged=converged,niter=fit_cox$iter,
               flag_nan=fit_cox$`NaN flag`)
    if(test){
      list(mets=mets,fit_cox=fit_cox)
    }else{
      mets
    }

    
    
  }#end foreach
  
  if(parallel){
    stopCluster(cl)
  }
  
  
  return(metrics)
  
}