#' @export
init = function(X, M=NULL, y, delta, k, alpha, lambda, eta, lambdaW, lambdaH,
                seed = 123,tol = 1e-6, imaxit = 30, verbose = TRUE, ninit = 20){

  X = as.matrix(X)
  if(is.null(M)){
    M = matrix(1,nrow=nrow(X),ncol=ncol(X))
  }

  p = nrow(X)
  n = ncol(X)
  c_best = 0
  for(i in 1:ninit){
    set.seed(seed*i)
    H0 = matrix(runif(n*k,0,max(X)),nrow=k)
    W0 = matrix(runif(p*k,0,max(X)),nrow=p)
    #beta0 = runif(k,-1,1)
    beta0 = rep(0,k)

    fit = optimize_loss_cpp(X, M, y, delta, W0, H0, beta0,
                            alpha, lambda, eta,
                            lambdaW, lambdaH,
                            tol, imaxit, verbose, TRUE)

    Z=t(X)%*%fit$W
    lp=Z%*%fit$beta
    c = cvwrapr::getCindex(lp,Surv(y,delta))
    if(c > c_best & !fit$nan_flag){#  & fit$convergence
      c_best=c
      fit_best=fit
    }
  }
  return(fit_best)
}
