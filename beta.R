#function which takes the regressor matrices Z and Yex, number of equations nEq
#number of observations p, mode of estimation "estMode" and number of iterations
#in the case of iterative feasible general least square

Beta <- function(mat.z
                 , mat.y.ex
                 , n.eq
                 , p
                 , est.mode
                 , iter){
  
  
  if(est.mode == "OLS"){                                                  #if OLS
    mat.beta <- solve(mat.z %*% t(mat.z), tol = 0) %*% mat.z %*% mat.y.ex     #solve the system and get the vector of betas
  }

  if(est.mode == "FGLS"){                                                  #if FGLS
  
    mat.beta <- solve(mat.z %*% t(mat.z), tol = 0) %*% mat.z %*% mat.y.ex     #solve the system and get the vector of betas
    mat.sigma <- Sigma(mat.z, mat.y.ex, mat.beta, n.eq)                   #get the covariance matrix of errors
 
    mat.omega <- kronecker(diag(p), mat.sigma)                            #get Omega
    mat.beta <- solve(mat.z %*% solve(mat.omega, tol = 0) %*% t(mat.z), tol = 0) %*% mat.z %*% solve(mat.omega, tol = 0) %*% mat.y.ex        #solve the system and get the vector of betas
  }

  if(est.mode == "IGLS"){                                                  #if IGLS

    mat.beta <- solve(mat.z %*% t(mat.z), tol = 0) %*% mat.z %*% mat.y.ex                      #solve the system and get the vector of betas
    
    for (i in 1:iter){
   
      mat.sigma <- Sigma(mat.z, mat.y.ex, mat.beta, n.eq)                           #get the covariance matrix of errors
      mat.omega <- kronecker(diag(p), mat.sigma)                                    #get Omega
      mat.beta <- solve(mat.z %*% solve(mat.omega, tol = 0) %*% t(mat.z), tol = 0) %*% mat.z %*% solve(mat.omega, tol = 0) %*% mat.y.ex      #solve the system and get the vector of betas
      # print(Sigma)                                                                #control for convergence / non divergence
    }
  }
  mat.sigma <- Sigma(mat.z, mat.y.ex, mat.beta, n.eq)                               #final estimation of Sigma
  return(list(mat.beta = mat.beta, mat.sigma = mat.sigma))
}
