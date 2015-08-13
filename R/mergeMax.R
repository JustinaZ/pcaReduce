#' Function for merging two clusters (max)
#' 
#' Takes in K number of clusters and outputs the indeces of clusters that will be merged. Clusters for merging are selected based on largest probability.
#' @name MergeM
#' @param dat  -- data given as matrix, where in rows are n cells, and in columns are q genes
#' @param cl_id  -- cluster indicator variable
#' @param cent  -- a matrix, where each row is centre/mean of each cluster
#' @param K  -- number of clusters
#' @return Returns indeces of clusters that should be merged
#' @export
#' 


MergeM <- function(dat, cl_id, cent, K){
  
  P <- mat.or.vec(K,K) # log values
  a <- combn(seq(1:K),2)
  v <- c()
  for (l in 1:(0.5*(K*K-K))){
    ind1 <- a[1,l]
    ind2 <- a[2,l]
    x1 <- dat[which(cl_id == ind1), ]
    x2 <- dat[which(cl_id == ind2), ]
    if (class(x1) =="numeric"){
      C1 <- diag(K-1)
      n1 <- 1
      mean1 <- x1
    }else{
      C1 <- cov(x1)
      n1 <- nrow(x1)
      mean1 <- colMeans(x1)
    }
    if (class(x2) =="numeric"){
      C2 <- diag(K-1)
      n2 <- 1
      mean2 <- x2
    }else{
      C2 <- cov(x2)
      n2 <- nrow(x2)
      mean2 <- colMeans(x2)
    }
    m1 <- cent[ind1,]
    m2 <- cent[ind2,]
    pi1 <- n1/(n1+n2)
    pi2 <- n2/(n1+n2)
    p2 <- DMnorm_max(rbind(x1,x2), (pi1*mean1 + pi2*mean2), (pi1*C1+pi2*C2), K )
    P[ind1, ind2] <- p2
    P[ind2, ind1] <- NA
    v <- c(v, p2)
  }
  diag(P) <- NA
  max_value <- which(P==max(v), TRUE)
  
  if (length(which(P==max(v), TRUE))/2 >1){
    to_merge <- max_value[sample(1:nrow(max_value),1), ]
  }else{
    to_merge <- max_value
  }
  
  c(to_merge, pi1, pi2)
}
