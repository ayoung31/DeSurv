#' @export
set_param_grid = function(k, alpha, lambda, eta, lambdaW, lambdaH,
                          ninit, type, prefix,
                          ngene, maxit, tol, imaxit,
                          nfold=NULL, nmask=NULL, perc_mask=NULL){
  
  if(type=="full"){
    if(!dir.exists(paste0("results/",prefix,"/full/ngene",ngene,"/raw/"))){
      dir.create(paste0("results/",prefix,"/full/ngene",ngene,"/raw/"),recursive = TRUE)
    }
    params = expand.grid(k=k,alpha=alpha,lambda=lambda,eta=eta,lambdaW=lambdaW,
                         lambdaH=lambdaH)

    params$file=paste0('results/',prefix,'/full/ngene',ngene,'/raw/k=',params$k,
                       '_alpha',params$alpha,'_lambda',params$lambda,'_eta',params$eta,
                       '_lambdaW',params$lambdaW, '_lambdaH',params$lambdaH,
                       '_full','_ninit',ninit,
                       '_imaxit',imaxit,'_tol',tol,'_maxit',maxit,'.RData')
  }else if(type=="cv"){
    if(!dir.exists(paste0("results/",prefix,"/cv/ngene",ngene,"/raw/"))){
      dir.create(paste0("results/",prefix,"/cv/ngene",ngene,"/raw/"),recursive = TRUE)
    }
    params = expand.grid(k=k,alpha=alpha,lambda=lambda,eta=eta,fold=1:nfold,
                         lambdaW=lambdaW,lambdaH=lambdaH)
    params$file=paste0('results/',prefix,'/cv/ngene',ngene,'/raw/k=',params$k,
                       '_alpha',params$alpha,'_lambda',params$lambda,'_eta',params$eta,
                       '_lambdaW',params$lambdaW, '_lambdaH',params$lambdaH,
                       '_fold',params$fold,'of',nfold,'_ninit',ninit,
                       '_imaxit',imaxit,'_tol',tol,'_maxit',maxit,'.RData')
  }else if(type=="mask"){
    if(!dir.exists(paste0("results/",prefix,"/mask/ngene",ngene,"/raw/"))){
      dir.create(paste0("results/",prefix,"/mask/ngene",ngene,"/raw/"),recursive = TRUE)
    }
    params = expand.grid(k=k,alpha=alpha,lambda=lambda,eta=eta,mask=1:nmask,
                         lambdaW=lambdaW,lambdaH=lambdaH)
    params$file=paste0('results/',prefix,'/mask/ngene',ngene,'/raw/k=',params$k,
                       '_alpha',params$alpha,'_lambda',params$lambda,'_eta',params$eta,
                       '_lambdaW',params$lambdaW, '_lambdaH',params$lambdaH,
                       '_percmask',perc_mask,'_mask',params$mask,'_ninit',ninit,
                       '_imaxit',imaxit,'_tol',tol,'_maxit',maxit,'.RData')
  }
  

  return(params)
}