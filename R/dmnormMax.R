#' Joint multivariate normal probability (max)
#' 
#' Takes in n data points and outputs their joint probability under multivariate normal distribution. This function is used by MergeM function.
#' @name DMnorm_max
#' @param x -- data given in a matrix format, where rows are n cells, and columns are q genes
#' @param m -- mean vector
#' @param C -- covariance matrix q x q
#' @param K -- number of clusters
#' @return Returns a joint probability of data (i.e. sum of log-density values)
#' 

DMnorm_max <- function(x, m, C, K){
  sum(dmnorm(x, m, (C+ 0.05*diag(K-1)), log=TRUE))
}

