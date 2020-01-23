compute_ci <- function(G, S, Sigma, R, Beta, cv, p)             # as per Bekeart Harvey Lumsdaine (2002) // recall that RBeta = S %*% deltaT
{
  u <- length(cv)                                               #get the number of critical values from the vdistr
  ciDelta <- rep(NA, u)                                         #create empty vector

  RBeta <- R %*% Beta                                           #pre compute the selection of (breaking) parameters to test

  tCi <- solve(t(RBeta) %*% S %*% kronecker((G %*% t(G)) / p, solve(Sigma, tol = 0)) %*% t(S) %*% RBeta, tol = 0)    #compute the conf interval factor

  for(i in 1:u)
  {
    ciDelta[i] <- cv[i] * tCi                                      #get the vector of critical values
  }
  # print(ciDelta)
  return(ciDelta)
}
