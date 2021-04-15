# as per Bekeart Harvey Lumsdaine (2002) // recall that mat.r.beta = mat.s %*% delta.t

ConfidenceInterval <- function(mat.g
                               , mat.s
                               , mat.sigma
                               , mat.r
                               , mat.beta
                               , cv
                               , p){             

  n <- length(cv)                                                                            #get the number of critical values from the vdistr
  ci.delta <- rep(NA, n)                                                                     #create empty vector

  mat.r.beta <- mat.r %*% mat.beta                                                           #pre compute the selection of (breaking) parameters to test

  t.ci <- solve(t(mat.r.beta) %*% mat.s %*% kronecker((mat.g %*% t(mat.g)) / p, solve(mat.sigma, tol = 0)) %*% t(mat.s) %*% mat.r.beta, tol = 0)    #compute the conf interval factor

  for(i in 1:n){
    
    ci.delta[i] <- cv[i] * t.ci                                                               #get the vector of critical values
  }
  return(ci.delta)
}
