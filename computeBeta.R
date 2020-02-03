compute_beta <- function(Z, Yex, nEq, p, estMode, iter)                   #function which takes the regressor matrices Z and Yex, number of equations nEq
  #number of observations p, mode of estimation "estMode" and number of iterations
  #in the case of iterative feasible general least square
{
  if(estMode == "OLS")                                                    #if OLS
  {
    Beta <- solve(Z %*% t(Z), tol = 0) %*% Z %*% Yex                      #solve the system and get the vector of betas
  }

  if(estMode == "FGLS")                                                   #if FGLS
  {
    Beta <- solve(Z %*% t(Z), tol = 0) %*% Z %*% Yex                      #solve the system and get the vector of betas
    Sigma <- compute_sigma(Z, Yex, Beta, nEq)                             #get the covariance matrix of errors
 
    Omega <- kronecker(diag(p), Sigma)                                    #get Omega
    Beta <- solve(Z %*% solve(Omega, tol = 0) %*% t(Z), tol = 0) %*% Z %*% solve(Omega, tol = 0) %*% Yex        #solve the system and get the vector of betas
  }

  if(estMode == "IGLS")                                                   #if IGLS
  {
    Beta <- solve(Z %*% t(Z), tol = 0) %*% Z %*% Yex                      #solve the system and get the vector of betas
    
    for (i in 1:iter)
    {
      Sigma <- compute_sigma(Z, Yex, Beta, nEq)                           #get the covariance matrix of errors
      Omega <- kronecker(diag(p), Sigma)                                  #get Omega
      Beta <- solve(Z %*% solve(Omega, tol = 0) %*% t(Z), tol = 0) %*% Z %*% solve(Omega, tol = 0) %*% Yex      #solve the system and get the vector of betas
      # print(Sigma)                                                      #control for convergence / non divergence
    }
  }
  Sigma <- compute_sigma(Z, Yex, Beta, nEq)                               #final estimation of Sigma
  return(list(Beta = Beta, Sigma = Sigma))
}
