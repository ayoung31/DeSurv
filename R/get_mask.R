#' @export
get_mask = function(perc_miss,n,p){
  M <- matrix(1,nrow=p,ncol=n)
  n_miss <- floor(n*p*perc_miss)
  entries <- as.vector(outer(1:n, 1:p, paste, sep="_"))
  miss_entries <- sample(entries,n_miss,replace=FALSE)
  rcs <- strsplit(miss_entries,'_')
  cols <- as.numeric(lapply(rcs,'[',1))
  rows <- as.numeric(lapply(rcs,'[',2))
  M[cbind(rows,cols)] <- 0
  
  return(M)
}
