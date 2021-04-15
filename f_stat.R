#compute the f-statistic as in BLS (1998)

Fstat <- function(mat.r
                  , mat.beta
                  , mat.z
                  , p
                  , mat.sigma){
  
  mat.r.beta <- mat.r %*% mat.beta                         #pre compute the RBeta matrix with the selected coefficients allowed to break

  if(!is.null(mat.sigma)){
    #if the covariance matrix of error is passed as argument, compute F-stat

    mat.omega <- kronecker(diag(p), mat.sigma)        #get Omega
    f.k <- p * t(mat.r.beta) %*% solve(mat.r %*% solve((mat.z %*% solve(mat.omega, tol = 0) %*% t(mat.z)) / p, tol = 0) %*% t(mat.r), tol = 0) %*% mat.r.beta
  }
  return(f.k)
}
