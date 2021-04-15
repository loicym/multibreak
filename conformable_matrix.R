

ConformableMatrix <- function(mat.y
                              , q
                              , mat.x
                              , trend
                              , intercept){
  
  p.init <-  dim(mat.y)[1]                                                           #get the original number of observations

  l.y <- Lags(mat.y = mat.y
              , q = q)                                                      #get the list of contemporaneous and lags objects

  mat.y <- l.y$mat.y                                                                     #get the matching dependent matrix
  mat.y.lag <- l.y$mat.y.lag                                                               #get the lagged dependent variables matrix
  my.dates <- l.y$my.dates                                                         #get the original matching rowname vector (of dates)
  
  p <- dim(mat.y)[1]                                                                #length of the matrix
  n.eq <- dim(mat.y)[2]                                                              #number of equations of the vAR

  print(p)
  id.n <- diag(n.eq)                                                             #identity matrix of the number of equations of the VAR

  mat.g.t <- as.matrix(cbind(rep(1, p), mat.y.lag))                                       #create a unique matrix transpose of G with one intercept and autoregressive terms
  
  n <- dim(mat.g.t)[2]                                                               #incremental number of regressors

  if(!is.null(mat.x)){                                                               #if additional covariates matrix is passed as argument
    
    if(p.init == dim(mat.x)[1]){                                                        #and if its size is equal to the original matrix Y
      
      mat.g.t <- cbind(mat.g.t, data.matrix(mat.x[(q + 1) : p.init, ]))                        #increment Gt by the contemporaneous covariate matrix X
      n <- dim(mat.g.t)[2]                                                           #increment the total number of regressors
    }
    
    else
      print("The number of observations of X does not match the one of Y")
  }

  if(trend){                                                                     #check if we add a trend
    
    mat.g.t <- cbind(mat.g.t, seq(1, p, by = 1))                                          #add trend
    n <- dim(mat.g.t)[2]                                                             #increment the total number of regressors
  }

  if(intercept){                                                                 #if only the intercept is allowed to break
    
    s <- t(data.matrix(c(rep(0, n))))                                           #create the selection vector
    s[1] <- 1                                                                   #vector only select the first element (intercept) to test the shift
    mat.s <- kronecker(s, id.n)                                                       #create the selection matrix
  }

  if(!intercept){                                                                #full parameters structural change estimation
                                                                                #get the full dimension of the test
    r <- n.eq * n
    mat.s <- diag(r)                                                                #identity matrix of the size of the test is the selection matrix
  }

  mat.g <- t(mat.g.t)                                                                    #transpose Gt to get G as in BLS

  mat.y.ex <- data.matrix(c(t(mat.y)))                                                   #Expend and vectorize Y
  
  mat.g.ex <- kronecker(t(mat.g), id.n)                                                    #Expend G
  return(list(mat.y.ex = mat.y.ex
              , mat.g.ex = mat.g.ex
              , p = p
              , mat.g = mat.g
              , mat.s = mat.s
              , my.dates = my.dates
              , mat.y = mat.y
              , n.eq = n.eq))
}
