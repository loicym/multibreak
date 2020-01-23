compute_fstat <- function(R, Beta, Z, p, Sigma = NULL)                   #compute the f-statistic as in BLS (1998)
{
  RBeta <- R %*% Beta                         #pre compute the RBeta matrix with the selected coefficients allowed to break

  if(is.null(Sigma))                          #if no previously covariance matrix of error is passed as argument, compute F-stat
  {
    Fk <- p * t(RBeta) %*% solve(R %*% solve( (Z %*% t(Z) / p), tol = 0) %*% t(R), tol = 0) %*% RBeta
  }

  if(!is.null(Sigma))                         #if the covariance matrix of error is passed as argument, compute F-stat
  {
    Omega <- kronecker(Sigma, diag(p))
    Fk <- p * t(RBeta) %*% solve(R %*% solve((Z %*% solve(Omega, tol = 0) %*% t(Z) / p), tol = 0) %*% t(R), tol = 0) %*% RBeta
  }
  return(Fk)
}
