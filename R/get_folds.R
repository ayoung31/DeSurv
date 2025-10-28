#' @export
get_folds = function(X,y,nfold,ngene,qnorm=FALSE){
  
  folds=caret::createFolds(y,k=nfold,list=FALSE)
  
  Xtrain = list()
  Xtest = list()
  for(i in 1:nfold){
    
    #train_data = gene_filter(X[,folds!=i],ngene)
    train_data = X[,folds!=i]
    
    # reference_quantiles=rowMeans(apply(train_data,2,sort))
    # train_data_normalized = normalize.quantiles(train_data,keep.names = TRUE)
    # subset each testing fold to genes kept in training fold
    # test_data = X[rownames(X) %in% rownames(train_data),folds==i]
    test_data = X[,folds==i]
    # test_data_normalized = normalize.quantiles.use.target(test_data,target = reference_quantiles)
    
    if(qnorm){
      # Xtrain[[i]] = train_data_normalized
      # Xtest[[i]] = test_data_normalized
    }else{
      Xtrain[[i]] = as.matrix(train_data)
      Xtest[[i]] = as.matrix(test_data)
    }
  }#end for loop
  
  return(list(Xtrain=Xtrain,Xtest=Xtest,folds=folds))
}
