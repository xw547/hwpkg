################################################################################
###########################     Algorithm Leveraging     #######################
################################################################################
#' Algorithm leveraging function
#'
#' A function that applies the method mentioned by (Ma, 2014) to do a leveraging
#' through two sampling methods, the uniform sampling and weighted sampling
#'
#'@param y The constant vector
#'@param X The regressor
#'@param size Sample size of leveraging sampling.
#'@param draws Number of draws of subsample
#'@param method The method to implement, can be weighted leveraging and uniform
#'leveraging. Default is "weighted". Details of the difference of the
#'
#'@return The leverage solution to the linear regressioin problem.
#'
#'@author Xiaohan Wang
#'
#'@export


algo_leverage <- function(
              y,
              X,
              size = NA,
              draws = NA,
              method ="weighted"
                          ){
  n = length(y)
### Parameter Check
  if (is.matrix(X)){
    if (length(y) != dim(X)[1]){
      warning("Dimension of input doesn't match")

    }
  }else{
    if (length(y) != length(X)){
      warning("Dimension of input doesn't match")

    }
  }


  if (size%%1 != 0 || size<0 || is.nan(size) ){
    warning("size should be a positive interger")

  }

  if (draws%%1 != 0 || draws<0 || is.nan(draws) ){
    warning("draws should be a positive interger")

  }

###Existing Functions
  package_weighted <- function(r,Pi_weighted,draws) {
    Beta_Hat_weighted_500 = c()
    for (i in 1:draws){
      index_2=sample(c(1:n), size = r, replace = T,prob=Pi_weighted)
      X_star_2=X[index_2]
      y_star_2=y[index_2]
      Phi_2=1/ (Pi_weighted[index_2]^{1/2})
      model_2 <- lm(y_star_2 ~ 0 + X_star_2 , weights=Phi_2)

      beta_hat_2=model_2$coefficients
      Beta_Hat_weighted_500[i]=beta_hat_2
    }
    return(Beta_Hat_weighted_500)
  }

  package_unif <- function(r,Pi_unif, draws) {
    Beta_Hat_unif_500 = c()
    for (i in 1:draws){
      index_1=sample(c(1:n), size = r, replace = T,prob=Pi_unif)
      X_star_1=X[index_1]
      y_star_1=y[index_1]
      Phi_1=1/ (Pi_unif[index_1]^{1/2})
      model_1 <- lm(y_star_1 ~ 0 + X_star_1 , weights=Phi_1)
      beta_hat_1=model_1$coefficients
      Beta_Hat_unif_500[i]=beta_hat_1
    }
    return(Beta_Hat_unif_500)
  }

  if (method == "uniform"){
    Pi = rep(1/n, n)
    return(package_unif(size, Pi,draws))
  } else if(method  == "weighted"){
    XTXinv = solve(t(X)%*%X)
    H =X%*%XTXinv%*%t(X)
    H_ii = diag(H)/sum(diag(H))
    Pi = H_ii
    return(package_weighted(size, Pi,draws))
  } else {
    warning("invalid method parameter")

  }
}

