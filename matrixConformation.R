matrix_conformation <- function(Y, q, X, trend, intercept)
{
  pInit <-  dim(Y)[1]                                                           #get the original number of observations

  lY <- compute_lags(Y, q)                                                      #get the list of contemporaneous and lags objects

  Y <- lY$Y                                                                     #get the matching dependent matrix
  YLAG <- lY$YLAG                                                               #get the lagged dependent variables matrix
  myDates <- lY$myDates                                                         #get the original matching rowname vector (of dates)
  
  p <- dim(Y)[1]                                                                #length of the matrix
  nEq <- dim(Y)[2]                                                              #number of equations of the vAR

  print(p)
  In <- diag(nEq)                                                               #identity matrix of the number of equations of the VAR

  Gt <- as.matrix(cbind(rep(1, p), YLAG))                                       #create a unique matrix transpose of G with one intercept and autoregressive terms
  
  n <- dim(Gt)[2]                                                               #incremental number of regressors

  if(!is.null(X))                                                               #if additional covariates matrix is passed as argument
  {
    if(pInit==dim(X)[1])                                                        #and if its size is equal to the original matrix Y
    {
      Gt <- cbind(Gt, data.matrix(X[(q + 1) : pInit, ]))                        #increment Gt by the contemporaneous covariate matrix X
      n <- dim(Gt)[2]                                                           #increment the total number of regressors
    }
    else
      print("The number of observations of X does not match the one of Y")
  }

  if(trend)                                                                     #check if we add a trend
  {
    Gt <- cbind(Gt, seq(1, p, by = 1))                                          #add trend
    n <- dim(Gt)[2]                                                             #increment the total number of regressors
  }

  if(intercept)                                                                 #if only the intercept is allowed to break
  {
    s <- t(data.matrix(c(rep(0, n))))                                           #create the selection vector
    s[1] <- 1                                                                   #vector only select the first element (intercept) to test the shift
    S <- kronecker(s, In)                                                       #create the selection matrix
  }

  if(!intercept)                                                                #full parameters structural change estimation
  {                                                                             #get the full dimension of the test
    r <- nEq * n
    S <- diag(r)                                                                #identity matrix of the size of the test is the selection matrix
  }

  G <- t(Gt)                                                                    #transpose Gt to get G as in BLS

  Yex <- data.matrix(c(t(Y)))                                                   #Expend and vectorize Y
  
  Gex <- kronecker(t(G), In)                                                    #Expend G
  return(list(Yex = Yex, Gex = Gex, p = p, G = G, S = S, myDates = myDates, Y = Y, nEq = nEq))
}
