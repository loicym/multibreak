compute_fstat <- function(R, Beta, Z, p, Sigma)                   #compute the f-statistic as in BLS (1998)
{
  RBeta <- R %*% Beta                         #pre compute the RBeta matrix with the selected coefficients allowed to break

  if(!is.null(Sigma))                         #if the covariance matrix of error is passed as argument, compute F-stat
  {
    Omega <- kronecker(diag(p), Sigma)        #get Omega
    Fk <- p * t(RBeta) %*% solve(R %*% solve((Z %*% solve(Omega, tol = 0) %*% t(Z)) / p, tol = 0) %*% t(R), tol = 0) %*% RBeta
  }
  return(Fk)
}
