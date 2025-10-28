#' @export
gene_filter = function(X, ngene=5000){

  means = apply(X,1,mean) # average expression for each genes
  Xtrain_temp = X[means > quantile(means,.25),] # take top 75% of genes based on mean expression
  stds = apply(Xtrain_temp,1,sd) # Take std of expression for each genes
  train_data = Xtrain_temp[stds>quantile(stds,1-(ngene/length(stds))),] # Keep top 5000 most variable genes
  
  return(train_data)
}