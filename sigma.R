#compute the covariance matrix of errors as in BLS (1998)

Sigma <- function(mat.z
                  , mat.y.ex
                  , mat.beta
                  , n.eq){                  

  errors <- (mat.y.ex - t(mat.z) %*% mat.beta)                               #get the n*p vector of residuals
  mat.errors <- matrix(errors, ncol = n.eq, byrow = T)                       #reset as matrix to obtain Sigma
  mat.sigma <- cov(mat.errors)

  return(mat.sigma)
}
