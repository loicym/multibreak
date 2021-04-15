#function returning a list of matching matrices of dependent (Y) and independent (YLAG) lagged observations

Lags <- function(mat.y = mat.y,
                 q = q){
  
  p <- dim(mat.y)[1]                                                              #get the dimensions
  n <- dim(mat.y)[2]

  my.dates <- rownames(mat.y)[(q + 1):p]                                          #optional: keep the rownames dates of the data frame with final matching
  mat.y <- data.matrix(mat.y)                                                          #matrix conversion

  mat.y.lag <- matrix(data = NA, nrow <- (p - q), ncol <- (n * (q + 1)))            #create an empty matrix

  for(i in 0:q){
    
    mat.y.lag[ , (n * i + 1):(n * (i + 1))] <- mat.y[(q - i + 1):(p - i), ]
  }

  mat.y <- data.matrix(mat.y.lag[, 1:n])
  mat.y.lag <- data.matrix(mat.y.lag[, (n + 1):dim(mat.y.lag)[2]])
  return(list(mat.y = mat.y, mat.y.lag = mat.y.lag, my.dates = my.dates))
}
