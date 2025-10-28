#' @export
recommend_alpha <- function(X,M,y,delta,k,nalpha,maxit=15,WtX,eta=0,lambda=0, norm.type){
  if(nalpha < 1){
    warning('nalpha must be a positive integer')
  }
  p <- nrow(X)
  n <- ncol(X)
  set.seed(12)
  H0 = matrix(runif(n*k,0,10),nrow=k)
  W0 = matrix(runif(p*k,0,1),nrow=p)
  W0 = W0 %*% diag(1/colSums(W0))
  beta0 = runif(k,-1,1)
  fit0 <- optimize_loss_cpp(X,M,H0,W0,beta0,y,delta,0,lambda,eta,1e-6,maxit,
                            FALSE,WtX,norm.type,'lasso',TRUE)
  
  print(fit0$loss$nmf_loss/sum(M))
  print(2*fit0$loss$surv_loss/n)
  if(fit0$loss$nmf_loss > fit0$loss$surv_loss){
    alpha5050 <- -1 * fit0$loss$nmf_loss / fit0$loss$surv_loss
  }else{
    alpha5050 <- -1 * fit0$loss$surv_loss / fit0$loss$nmf_loss
  }
  if(nalpha %% 2 > 0){
    nalpha=nalpha+1
  }
  grid <- 2:(nalpha/2 + 1)
  alpha_grid <- c(alpha5050*grid,alpha5050,alpha5050/grid)
  alpha_grid = alpha_grid[order(alpha_grid)]
  alpha_grid = round(alpha_grid,3)
  alpha_grid = alpha_grid[order(alpha_grid)]
  return(list(alpha5050=alpha5050, alpha_grid=alpha_grid))
}
