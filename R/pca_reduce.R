#'Function to perform clustering and dimensionality reduction
#'
#'Takes in your dataset (a matrix with rows corresponding to cells, and columns corresponding to genes), performs prespecified dimensionality reduction and merning.
#'@name PCAreduce
#'@param D_t -- dataset, presented as marix, where rows are cells, and colms are genes
#'@param nbt -- number of samples, i.e. number of times to repeat pcaReduce framework
#'@param q -- number of dimensions to start with
#'@param method -- can be S or M, S - to perform sampling based merging, M - to perform merging based on largest probability
#'@return Returns a list of length nbt; each item in the list is a matrix containing cell allocation variables 
#'@export

PCAreduce <- function(D_t, nbt, q, method){
  Y <- prep(D_t, scale="none", center=TRUE)
  pca_out <- pca(Y, method="svd", center = FALSE, nPcs=q)
  x <- pca_out@scores
  
  if (method == "S"){
    Cl_history <- list() # nbt length
    for (t in 1:nbt){
      dat <- x
      q <- ncol(dat)
      K <- q + 1
      KM <- kmeans(dat, K)
      cent <- KM$centers
      cl_id <- KM$cluster
      
      Cl_mat <- c(cl_id)
      Q <- q-1
      for (i in 1:Q){
        mrg <- MergeS(dat, cl_id, cent, K)
        cl_id[which(cl_id==mrg[[1]][2])] <- mrg[[1]][1]
        cent[mrg[[1]][1],] <- mrg[[2]][1] * cent[mrg[[1]][1],] + mrg[[2]][2] * cent[mrg[[1]][2],]
        cent <- cent[-mrg[[1]][2],-ncol(cent)]
        
        a <- unique(cl_id)
        Omega <- seq(1, max(a), 1)
        b <- setdiff(Omega, a)
        N <- length(b)
        if (N>0){
          for (ii in 1:N){
            cl_id[cl_id>b[N+1-ii]] <- cl_id[cl_id>b[N+1-ii]]-1
          }
        }
        K <- length(unique(cl_id))
        dat <- dat[,-ncol(dat)]
        
        Cl_mat <- cbind(Cl_mat, cl_id)
      }
      Cl_history[[t]] <- Cl_mat
    }
    return(Cl_history)
  }else if (method == "M"){
    Cl_history<- list()
    for (t in 1:nbt){
      dat <- x
      q <- ncol(dat)
      K <- q + 1
      KM <- kmeans(dat, K)
      cent <- KM$centers
      cl_id <- KM$cluster
      Cl_mat <- c(cl_id)
      Q <- q - 1
      for (i in 1:Q){
        mrg_ind <- MergeM(dat, cl_id, cent, K)
        #cat("merge=", mrg_ind,"\n")
        cl_id[which(cl_id==mrg_ind[2])] <- mrg_ind[1]
        cent[mrg_ind[1],] <- mrg_ind[3] * cent[mrg_ind[1],] + mrg_ind[4] * cent[mrg_ind[2],]
        cent <- cent[-mrg_ind[2],-ncol(cent)]
        
        a <- unique(cl_id)
        Omega <- seq(1, max(a), 1)
        b <- setdiff(Omega, a)
        N <- length(b)
        if (N>0){
          for (ii in 1:N){
            cl_id[cl_id>b[N+1-ii]] <- cl_id[cl_id>b[N+1-ii]]-1
          }
        }
        K <- length(unique(cl_id))
        dat <- dat[,-ncol(dat)]
        Cl_mat <- cbind(Cl_mat, cl_id)
      }
      Cl_history[[t]] <- Cl_mat
    }
    return(Cl_history)
  }
  
}