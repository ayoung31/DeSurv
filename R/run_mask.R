#' @export
run_mask = function(X, y, delta, k, alpha, lambda=0, eta=0, 
                    lambdaW=0, lambdaH=0,
                    perc_mask, nmask, 
                    ninit=100, imaxit=30, maxit=3000, tol=1e-5,
                    parallel=TRUE,ncore=NULL,replace=FALSE,
                    save=TRUE, verbose=TRUE, prefix){
  
  X=as.matrix(X)
  
  if(ncore==1){
    parallel=FALSE
  }
  
  if(parallel & is.null(ncore)){
    ncore = detectCores() - 1
  }
  
  params = set_param_grid(k=k, alpha=alpha, lambda=lambda, eta=eta, ninit=ninit,
                          lambdaW = lambdaW, lambdaH = lambdaH,
                          type="mask", perc_mask=perc_mask, nmask=nmask, 
                          prefix=prefix, ngene = ngene, maxit=maxit, tol=tol, imaxit=imaxit)
  
  n=ncol(X)
  p=nrow(X)
  
  if(parallel){
    cl = parallel::makeCluster(ncore,outfile="")
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  }

  metrics = foreach(pa=1:nrow(params), .inorder = FALSE, .errorhandling = 'remove', 
                    .combine = 'rbind', .packages=c('coxNMF',"survival","cvwrapr")) %dopar% {
    
    
    a = params$alpha[pa]
    l = params$lambda[pa]
    e = params$eta[pa]
    k = params$k[pa]
    m = params$mask[pa]
    lW = params$lambdaW[pa]
    lH = params$lambdaH[pa]
    
    if(replace | !file.exists(params$file[pa])){
      set.seed(m)
      M = get_mask(perc_mask,n,p)
      print(sprintf('k=%d, alpha=%.2f, lambda=%.4f, eta=%.2f, lambdaW=%.1f, lambdaH=%.5f, mask=%d\n',k,a,l,e,lW,lH,m))
      print("running...")
      
      fit_cox = run_coxNMF(X=X, y=y, delta=delta, k=k, M=M,
                           alpha=a, lambda=l, eta=e, 
                           lambdaW = lW, lambdaH = lH,
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
    

    converged=fit_cox$iter < maxit
    
    W = fit_cox$W
    H = fit_cox$H
    beta = fit_cox$beta
    
    if(any(is.nan(t(X) %*% W %*% beta))){
      warning("alpha too large 2")
    }
    
    c = cvwrapr::getCindex(t(M*X) %*% W %*% beta, Surv(y, delta))
    cm = cvwrapr::getCindex(t((1-M)*X) %*% W %*% beta, Surv(y, delta))
    loss = fit_cox$loss
    lossm = calc_loss_cpp(X,(1-M),y,delta,W,H,beta,a,l,e,1,1,lW,lH)
    ol = loss$loss
    olm = lossm$loss
    sl = loss$surv_loss
    slm = lossm$surv_loss
    nl = loss$nmf_loss
    nlm = lossm$nmf_loss
    bic = -2*sl + k*log(ncol(X))
    bicm = -2*slm + k*log(ncol(X))
    
    data.frame(k=k,alpha=a,lambda=l,eta=e,lambdaW=lW, lambdaH=lH, mask=m,
               c=c,cm=cm,loss=ol,lossm=olm,sloss=sl,slossm=slm,nloss=nl,nlossm=nlm,
               bic=bic,bicm=bicm,penW=loss$penalty_W,penH=loss$penalty_H,
               converged=converged,niter=fit_cox$iter,
               flag_nan=fit_cox$`NaN flag`)
    
  }#end foreach
  
  if(parallel){
    stopCluster(cl)
  }
  
  return(metrics)
  
}