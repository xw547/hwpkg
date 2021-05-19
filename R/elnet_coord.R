################################################################################
#########################          Elastic Net         #########################
################################################################################

#' Elastic Net
#'
#' A function that implements the algorithm proposed by (Zou, Hui; Hastie,
#' Trevor, 2015, JRRSB).
#'
#'@param y The constant vector
#'@param X The regressor
#'@param beta0 The initial guess of the solution. If no input given, then it is
#'set to a zero vector.
#'@param alpha The elastic net parameter, alpha should be in [0,1]. More details
#'to be found in Trevors paper.
#'@param niter Maximum Number of Iteration
#'@param tol Tolerance of the difference between two result from iteration
#'
#'@return The solution to the elastic net system according to given parameter.
#'
#'@author Xiaohan Wang
#'
#'@export



elnet_coord <- function( y,
                         X,
                         beta0 = NA,
                         alpha = 0,
                         niter = 1e+4,
                         tol = 1e-5
){
###Parameter Check
  if (is.matrix(X)){
    if (length(y) != dim(X)[1]){
      warning("Dimension of input doesn't match")

    }
  }else{
    if (length(y) != length(X)){
      warning("Dimension of input doesn't match")

    }
  }
  if (alpha<0 || alpha>1 || isnan(alpha)){
    warning("the value of alpha isn't valid, should be within [0,1]")

  }
  if (niter%%1 != 0 || niter<0 || is.nan(niter) ){
    warning("niter should be a positive interger")

  }

  if (tol<0 || is.nan(tol) ){
    warning("tol should be a positive number")

  }

  if(is.na(beta0)){
    beta0 = rep(1, dim(X)[2])
  }
  return(Elastic(y, X, beta0, alpha, lambda, niter, tol))
}


###Existing functions
#'@export
soft_threshhold <- function(innerprod, lambda) {
  return(sign(innerprod)*max(0, abs(innerprod)- lambda))
}

#'@export
update_beta <- function(y, X, beta_new, lambda, alpha){
  beta_new_2 =c()
  n = length(y)
  a = colSums(X^2)
  p = length(beta_new)
  for(j in 1:p){
    beta_new[j] = soft_threshhold((crossprod(X[,j], (y-X[,-j]%*%beta_new[-j])))/n, lambda*alpha)
    beta_new[j] = beta_new[j]/(a[j]/n+lambda*(1-alpha))
  }
  return(beta_new)
}

#'@export
Elastic <- function(y, X, beta0, alpha, lambda, niter = 1e+4, tol = 1e-5){
  error = 1000
  n = 1
  all.beta = beta0
  gamma = lambda*alpha
  a = colSums(X^2)
  while (error > tol) {
    beta_new <- c(update_beta(y, X, beta0, lambda, alpha))
    all.beta <- rbind(all.beta, beta_new)
    if (any(abs(beta_new) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(beta_new - beta0)^0.5
    if (n == niter) {
      warning("Max iteration reached")
      break
    }
    beta0 <- beta_new
    n <- n + 1
  }
  return(beta_new)
}






