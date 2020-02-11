compute_beta_one <- function(Z        #function which takes the regressor matrices Z and Yex, number of equations nEq
                             , Yex    #number of observations p, mode of estimation "estMode" and number of iterations 
                             , nEq    #in the case of iterative feasible general least square
                             , p
                             , estMode
                             , iter) 
                   
  
{
  Yo <- matrix(Yex, ncol = nEq, byrow = T)   #reselect equations one by one to compute Betas and covariance matrix
  Zo <- matrix(Z, nrow = p, byrow = T)

  lZ <- list()
  mySeqRow <- seq(1, dim(Z)[1], by = nEq)
  mySeqCol <- seq(1, dim(Z)[2], by = nEq)

  Zo <- Z[mySeqRow, mySeqCol]
  Betao <- matrix(data = NA, nrow = dim(Zo)[1], ncol = nEq)
  Eo <- matrix(data = NA, nrow = dim(Yo)[1], ncol = dim(Yo)[2])
  
  for(i in 1:nEq)
  {
    Betao[, i] <- solve(Zo %*% t(Zo), tol = 0) %*% Zo %*% Yo[, i]
    Eo[,i] <- Yo[, i] - t(Zo) %*% Betao[, i]
  }
  
  Betao = matrix(c(t(Betao)), ncol = 1)

  Sigmao <- cov(Eo)                                                       #final estimation of Sigma
  return(list(Beta = Betao, Sigma = Sigmao))
}
