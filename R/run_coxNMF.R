#' @export
run_coxNMF = function(X, y, delta, k, alpha, M = NULL, W0 = NULL, H0 = NULL,
                      beta0 = NULL, lambda = 0, eta = 0,
                      lambdaW = 0, lambdaH = 0, seed = 123,
                      tol = 1e-6, maxit = 3000, verbose = TRUE, 
                      ninit = 30, imaxit = 100){
  # Add all error checking here. This will be the primary function
  genes = rownames(X)
  X = as.matrix(X)
  
  if(is.null(M)){
    M = matrix(1,nrow=nrow(X),ncol=ncol(X))
  }
  
  # Initialize
  if(is.null(H0) | is.null(W0) | is.null(beta0)){
    fit = init(X=X, M=M, y=y, delta=delta, k=k, 
                alpha=alpha, lambda=lambda, eta=eta,
                lambdaW, lambdaH,seed = seed,
                verbose=verbose, tol=tol, imaxit=imaxit, ninit=ninit)
  }else{
    # Run the model
    fit = optimize_loss_cpp(X, M, y, delta, W0, H0, beta0, 
                            alpha, lambda, eta, lambdaW, lambdaH,
                            tol, maxit, verbose, FALSE)
  }
  fit$genes = genes
  rownames(fit$W) = genes

  
  return(fit)
}
