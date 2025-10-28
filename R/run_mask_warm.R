#' @export
run_mask_warm = function(X, y, delta, k, alpha, lambda=0, eta=0, 
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
  alpha=alpha[order(alpha)]
  if(alpha[1]!=0){
    alpha=c(0,alpha)
  }
  
  params = set_param_grid(k=k, alpha=alpha, lambda=lambda, eta=eta, ninit=ninit,
                          lambdaW = lambdaW, lambdaH = lambdaH,
                          type="mask", perc_mask=perc_mask, nmask=nmask, 
                          prefix=prefix, ngene = ngene, maxit=maxit, tol=tol, imaxit=imaxit)
  
  params2=params %>% filter(alpha==0)
  
  n=ncol(X)
  p=nrow(X)
  
  if(parallel){
    cl = parallel::makeCluster(ncore,outfile="")
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  }

  metrics = foreach(pa=1:nrow(params2), .inorder = FALSE, .errorhandling = 'remove', 
                    .combine = 'rbind', .packages=c('coxNMF',"survival","cvwrapr")) %dopar% {
    
    l = params2$lambda[pa]
    e = params2$eta[pa]
    k = params2$k[pa]
    m = params2$mask[pa]
    lW = params2$lambdaW[pa]
    lH = params2$lambdaH[pa]
    
    set.seed(m)
    M = get_mask(perc_mask,n,p)
    
    
    dat = list()
    i=1
    for(a in alpha){
      file = params$file[params$alpha==a & params$lambda==l & 
                           params$eta==e & params$k==k & params$lambdaW==lW &
                           params$lambdaH==lH & params$mask==m]
      print(sprintf('k=%d, alpha=%.2f, lambda=%.4f, eta=%.2f, lambdaW=%.1f, lambdaH=%.5f, mask=%d\n',k,a,l,e,lW,lH,m))
      if(replace | !file.exists(file)){
        
        print("running...")
        if(a==0){
          fit_cox = run_coxNMF(X=X, y=y, delta=delta, k=k, M=M,
                               alpha=a, lambda=l, eta=e, 
                               lambdaW = lW, lambdaH = lH,
                               tol=tol, maxit=maxit, verbose=verbose,
                               ninit=ninit, imaxit=imaxit)
        }else{
          fit_cox = run_coxNMF(X=X, y=y, delta=delta, k=k, M=M,
                               alpha=a, lambda=l, eta=e, 
                               lambdaW = lW, lambdaH = lH,
                               tol=tol, maxit=maxit, verbose=verbose,
                               ninit=ninit, imaxit=imaxit,
                               W0=fit_cox$W, H0=fit_cox$H, beta0=fit_cox$beta)
        }
        
        if(save){
          save(fit_cox,file=file)
        }
      }else{
        print("loading...")
        load(file)
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
      
      dat[[i]] = data.frame(k=k,alpha=a,lambda=l,eta=e,lambdaW=lW, lambdaH=lH, mask=m,
                 c=c,cm=cm,loss=ol,lossm=olm,sloss=sl,slossm=slm,nloss=nl,nlossm=nlm,
                 bic=bic,bicm=bicm,penW=loss$penalty_W,penH=loss$penalty_H,
                 converged=converged,niter=fit_cox$iter,
                 flag_nan=fit_cox$`NaN flag`)
      i = i+1
    }#end alpha loop
    
    do.call('rbind',dat)
    
  }#end foreach
  
  if(parallel){
    stopCluster(cl)
  }
  
  return(metrics)
  
}