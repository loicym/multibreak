compute_sigma <- function(Z, Yex, Beta, nEq)                    #compute the covariance matrix of errors as in BLS (1998)
{
  errors <- (Yex - t(Z) %*% Beta)                               #get the n*p vector of residuals
  Errors <- matrix(errors, ncol = nEq, byrow = T)               #reset as matrix to obtain Sigma
  Sigma <- cov(Errors)

  return(Sigma)
}
