compute_lags <- function(Y, q)                                                #function returning a list of matching matrices of dependent (Y) and independent (YLAG) lagged observations
{
  p <- dim(Y)[1]                                                              #get the dimensions
  n <- dim(Y)[2]

  myDates <- rownames(Y)[(q + 1) : p]                                          #optional: keep the rownames dates of the data frame with final matching
  Y <- data.matrix(Y)                                                          #matrix conversion

  YLAG <- matrix(data = NA, nrow <- (p - q), ncol <- (n * (q + 1)))            #create an empty matrix

  for(i in 0:q)
  {
    YLAG[ , (n * i + 1):(n * (i + 1))] <- Y[(q - i + 1):(p - i), ]
  }

  Y <- YLAG[,1:n]
  YLAG <- YLAG[,(n + 1):dim(YLAG)[2]]
  return(list(Y = Y, YLAG = YLAG, myDates = myDates))
}
