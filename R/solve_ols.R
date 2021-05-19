################################################################################
########################### Solve Ordinary Linear System #######################
################################################################################
#'Solve linear system
#'
#' solve.ols is a function that solves linear system using either Gauss-Seidel
#' or Jacobi. It assumes that the spectral norm of A is less than 1 to solve.
#'
#'@param A Linear System to be solved
#'@param b Constant Vector
#'@param x0 Initial solution guess
#'@param maxiter Maximum Number of Iteration
#'@param tol Tolerance of the difference between two result from iteration
#'@param parallel Whether to use the paralleled version, Only available when
#'using "Jacobi". Default is FALSE.
#' @param method Method to use, default is "Gauss-Seidel", and can be "Jacobi"
#' @param full_result Whether to return a result of each iteration
#' @param numCores Only applicable when parallel == TRUE, this
#'
#'@return The function would return the solution to the linear system or the
#'full record of each iteration. This depends on the option of full_result.
#'
#'@author Xiaohan Wang
#'
#'@export
solve_ols <- function(A,                #Linear System to be solved
                      b,                #Constant Vector
                      x0 = NA,          #Initial solution
                      maxiter = 10000,  #Maximum Number of Iteration
                      tol =1e-5,        #Tolerance
                      numCores = 1,     #Number of cores when using the parallel
                      parallel = FALSE, #Whether to use the paralled version,
                                        #Only available when using "para"
                                        #Default is FALSE
                      method = "Gauss", #Method to use, default is "Jacobi"
                      full_result = FALSE
                      ){
###Parameter Check
if (is.matrix(A)){
    if (length(b) != dim(A)[1]){
      warning("Dimension of input doesn't match")

    }
  }else{
    if (length(b) != length(A)){
      warning("Dimension of input doesn't match")

    }
  }


if (eigen(solve(diag(diag(A)), getutri(A)+getltri(A)))$values[1] >= 1){
  warning("The problem can't be solved with Jacobi or Gauss-Seidel")

}

if (dim(A)[1]!= length(b)){
  warning("The design and output matrix isn't of same dimension")

}


if (is.na(x0)){
x0 = rep(0, dim(A)[2])
}

if(method == "Gauss" && parallel){
  warning("The parallel is only applicable when using the Jacobi method")
}

###Solve the problem using designated method.
  if (method == "Gauss"){
    result = GSMethod(A, b, x0, maxiter, tol)
  }
  if (method == "Jacobi"){
    if(parallel){
    result = seq_Jacobi(A, b, x0, maxiter, tol)
    }else{
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    result = para_Jacobi(A, b, x0, maxiter, tol)
    }
  }
  if (full_result){
    return(result)
  }else{
    return(tail(result, n = 1))
  }
}

#'@export
###The helper and the actual functions to implement
getutri <- function(A){
  return(as.matrix(A-tril(A)))
}

#'@export
###Same idea from the previous one.
getltri <- function(A){
  return(as.matrix(A-triu(A)))
}

#'@export
###Gauss-Seidel Method
GS <- function(A, b, x) {
  a <- diag(A)
  diag(A) <- 0
  for (i in 1:length(x)) {
    x[i] <- (b[i] - crossprod(A[i, ], x))/a[i]
  }
  return(x)
}

#'@export
GSMethod <- function(A, b, x0, maxiter = 10000, tol =1e-5){
  error = 1000
  n = 1
  all.x = x0
  while (error > tol) {
    x <- c(GS(A, b, x0))
    all.x <- rbind(all.x, x)
    if (any(abs(x) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(x - x0)^0.5
    if (n == maxiter) {
      warning("Max iteration reached")
      break
    }
    n <- n + 1
    x0 <- x
  }
  return(all.x)
}


###Jacobi Method
#'@export
Jacobi <- function(A, b, x){
  a <- diag(A)
  diag(A) <- 0
  x_new <- c()
  for (i in 1:length(x)) {
    x_new[i] <- (b[i] - crossprod(A[i, ], x))/a[i]
  }
  return(x_new)
}

#'@export
seq_Jacobi <- function(A, b, x0, maxiter = 10000, tol =1e-5){
  error = 1000
  n = 1
  all.x = x0
  while (error > tol) {
    x <- c(Jacobi(A, b, x0))
    all.x <- rbind(all.x, x)
    if (any(abs(x) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(x - x0)^0.5
    if (n == maxiter) {
      warning("Max iteration reached")
      break
    }
    n <- n + 1
    x0 <- x
  }
  return(all.x)
}

#'@export
###Parallel Jacobi
PJacobi <- function(A, b, x){
  a <- diag(A)
  diag(A) <- 0
  outlist = foreach (i = 1:length(x),  .multicombine = TRUE)  %dopar% {
    (b[i] - crossprod(A[i, ], x))/a[i]
  }
  x_new = as.vector(unlist(outlist))
  return(x_new)
}
#'@export
para_Jacobi <- function(A, b, x0, maxiter = 10000, tol =1e-5){
  error = 1000
  n = 1
  all.x = x0
  while (error > tol) {
    x <- c(PJacobi(A, b, x0))
    all.x <- rbind(all.x, x)
    if (any(abs(x) == Inf))
      stop("The algorithm diverges")
    error <- crossprod(x - x0)^0.5
    if (n == maxiter) {
      warning("Max iteration reached")
      break
    }
    n <- n + 1
    x0 <- x
  }
  return(all.x)
}

