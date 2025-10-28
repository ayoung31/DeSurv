#' @export
run_cv_warm = function(X, y, delta, k, nfold, alpha, lambda = 0, eta = 0, 
                  lambdaW = 0, lambdaH = 0, fold_info, seed = 123,
                  ninit = 100, imaxit=5000, maxit = 5000, tol = 1e-6, 
                  parallel = TRUE, ncore = NULL, replace = FALSE, 
                  save = TRUE, verbose=TRUE, prefix){
  
  X = as.matrix(X)
  
  alpha = alpha[order(alpha,decreasing = FALSE)]
  
  if(alpha[1] != 0){
    alpha = c(0,alpha)
  }
  
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
  
  params2=params %>% filter(alpha==0)
  
  Xtrain = fold_info$Xtrain
  Xtest = fold_info$Xtest
  folds = fold_info$folds
  
  if(parallel){
    cl = parallel::makeCluster(ncore,outfile="")
    doParallel::registerDoParallel(cl)
    parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  }
  
  metrics = 
    foreach(pa=1:nrow(params2), 
            .inorder = FALSE, 
            .errorhandling = 'remove', 
            .combine = 'rbind',
            .packages = c("coxNMF","survival","cvwrapr")) %dopar% 
    {#begin foreach
    
    suppressWarnings(rm(fit_cox))
    
    
    #a = params$alpha[pa]
    l = params2$lambda[pa]
    e = params2$eta[pa]
    f = params2$fold[pa]
    k = params2$k[pa]
    lW = params2$lambdaW[pa]
    lH = params2$lambdaH[pa]
    
    Xtr = Xtrain[[f]]
    y_curr = y[folds!=f]
    delta_curr = delta[folds!=f]
    
    dat = list()
    j=1
    for(a in alpha){
      file = params$file[params$alpha==a & params$lambda==l & 
                           params$eta==e & params$k==k & params$lambdaW==lW &
                           params$lambdaH==lH & params$fold==f]
      
    
      if(replace | !file.exists(file)){
        print(sprintf('alpha=%f, lambda=%f, eta=%f, k=%d, f=%d\n',a,l,e,k,f))
        print("running...")
        
        if(a==alpha[1]){
          
          p = nrow(X)
          n = ncol(X)
          loss_best = Inf
          suppressWarnings(rm(fit_best))
          for(i in 1:ninit){
            set.seed(seed*i)
            H0 = matrix(runif(n*k,0,max(X)),nrow=k)
            W0 = matrix(runif(p*k,0,max(X)),nrow=p)
            #beta0 = runif(k,-1,1)
            beta0 = rep(0,k)
            
            #subset H0 for current fold
            H0 = H0[,fold_info$folds!=f]
            M = matrix(1,nrow=nrow(Xtr),ncol=ncol(Xtr))
            
            
            fit = optimize_loss_cpp(Xtr, M, y_curr, delta_curr, W0, H0, beta0, 
                                    a, l, e,
                                    lW, lH,
                                    tol, imaxit, verbose, TRUE)
            
            loss = -1*fit$loss$surv_loss
            if(loss < loss_best & !fit$`NaN flag`){#  & fit$convergence
              loss_best=loss
              fit_best=fit
            }
            print(i)
          }
          
          fit_cox = fit_best
        }else{
          fit_cox = run_coxNMF(X=Xtr, y=y_curr, delta=delta_curr, k=k,
                               alpha=a, lambda=l, eta=e,
                               lambdaW = lW, lambdaH = lH,seed = seed,
                               tol=tol, maxit=maxit, verbose=verbose,
                               ninit=ninit, imaxit=imaxit,
                               W0=fit_cox$W, H0=fit_cox$H, beta0=fit_cox$beta)
        }
        
        if(save){
          save(fit_cox,file=file)
        }
      }else{
        print("loading ...")
        load(file)
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
      if(fit_cox$`NaN flag`){
        converged=FALSE
      }
      
      dat[[j]] = data.frame(k=k,alpha=a,lambda=l,eta=e,lambdaW=lW, lambdaH=lH,
                 fold=f,sloss=sl,sall=slall,iter=fit_cox$iter,
                 strain=fit_cox$loss$surv_loss,bic=bic,c=c, cvl= cvl, 
                 converged=converged, flag_nan=fit_cox$`NaN flag`)
    
      j=j+1
    }#end alpha loop
    do.call('rbind',dat)
  }#end foreach
  
  if(parallel){
    stopCluster(cl)
  }

  return(metrics)

}
