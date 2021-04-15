#function which takes the regressor matrices Z and Yex, number of equations nEq
#number of observations p, mode of estimation "estMode" and number of iterations 
#in the case of iterative feasible general least square

BetaOne <- function(mat.z                                  
                     , mat.y.ex                    
                     , n.eq                        
                     , p
                     , est.mode
                     , iter) {
  
  mat.y0 <- matrix(mat.y.ex, ncol = n.eq, byrow = T)   #reselect equations one by one to compute Betas and covariance matrix

  l.z <- list()
  seq.row <- seq(1, dim(mat.z)[1], by = n.eq)
  seq.col <- seq(1, dim(mat.z)[2], by = n.eq)

  mat.z0 <- mat.z[seq.row, seq.col]
  mat.betao <- matrix(data = NA, nrow = dim(Z0)[1], ncol = n.eq)
  e0 <- matrix(data = NA, nrow = dim(mat.y0)[1], ncol = dim(mat.y0)[2])
  
  for(i in 1:n.eq){
    
    mat.beta0[, i] <- solve(mat.z0 %*% t(mat.z0), tol = 0) %*% mat.z0 %*% mat.y0[, i]
    e0[, i] <- mat.y0[, i] - t(mat.z0) %*% mat.beta0[, i]
  }
  
  mat.beta0 <- matrix(c(t(mat.beta0)), ncol = 1)

  mat.sigma0 <- cov(e0)                                                       #final estimation of Sigma
  return(list(mat.beta = mat.beta0, mat.sigma = mat.sigma0))
}
