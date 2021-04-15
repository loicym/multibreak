#main entry function to compute BLS (1998) test of structural break.

#Y is a matrix or vector which will be lagged by (q) to compute a VAR(q)
#X is a matrix of (contemporaneous) covariates
#trend is a boolean indicating whether a trend vector should be added to the VAR
#intercept is a boolean indicating whether the test applies on the mean shift (TRUE) or all parameters (FALSE)
#ci is the vector of confidence intervals (in growing order) to compute based on the CDF of a V distr.
#estMode can take values of "OLS", "FGLS", "IGLS"
#in the case of "IGLS", the number of iteration "iter" can be specified.
#AicbicMode can be "AIC" or "BIC" depending on the minimum criterion to select
#qMax is the number of lags (from 1 to qMax) tested to determine the AIC / BIC maximum.
#trim is the percentage parameter to start and end the sample analysis
#if we want the algorithm to only detect positive breaks

Main <- function(mat.y                                      
                 , mat.x = NULL                             
                 , trend = FALSE                        
                 , intercept = TRUE                     
                 , ci = c(0.9, 0.95, 0.99)              
                 , est.mode = "OLS"                      
                 , iter = 3                             
                 , aic.bic.mode = "AIC"                   
                 , q.max = 6                             
                 , trim = 0.15                          
                 , pos.break = FALSE){
  
  my.vars <- colnames(mat.y)                                                   #get the variable names

  q.opt <- AicBic(mat.y, q.max, mat.x, trend, intercept)                    #return the AIC and BIC criteria for lags from 1 to qMax
  
  
  print(q.opt)                                                             #print the matrix of AIC and BIC for each lags

  q <- as.numeric(which.min(q.opt[aic.bic.mode, ]))                        #choose the lag q according to the min AIC
  print(paste0("lag with the minimum ", aic.bic.mode, " is: ", q))

  l.conf <- ConformableMatrix(mat.y = mat.y
                              , q = q
                              , mat.x = mat.x
                              , trend = trend
                              , intercept = intercept)                    #create a list of conform objects for the estimation

  mat.y.ex <- l.conf$mat.y.ex                                                        #get the conformed (expanded) Yex matrix (for the system, in vector form)
  mat.g.ex <- l.conf$mat.g.ex                                                        #get the conformed G matrix of regressors for the system
  p <- l.conf$p                                                            #final number of observations
  mat.g <- l.conf$mat.g                                                            #original matrix of regressors
  mat.s <- l.conf$mat.s                                                            #selection matrix
  mat.y <- l.conf$mat.y                                                            #matching original dependent variables matrix
  n.eq <- l.conf$n.eq                                                        #original number of equations / dependent variables
  my.dates <- l.conf$my.dates                                                #matching dates

  print(paste0("The number of equations in the system is: ", n.eq))

  f.stat <- rep(NA, p)                                                     #create a vector of f_statistics for each k tested
  mean.shift <- rep(NA, p)                                                 #create a vector with the evaluated size of the intercept difference
  mat.ci <- matrix(data = NA, nrow = p, ncol = length(ci))                 #create a matrix of confidence intervals for each k tested

  cv <- Vdistr(ci)                                                         #compute critical values for the vector of confidence intervals proposed

  start.ind <- round(trim * p)                                             #start index
  end.ind <- round(p - trim * p)                                           #end index

  for(k in start.ind:end.ind){                                             #loop over the k with a trimming date / burn period
  
    if(k%%10==0)
      print(paste0("The iteration is at the level k = ", k))            #get an idea of where we are in the loop every 10 iterations
    
    
    mat.g.ex.b <- mat.g.ex %*% t(mat.s)
    mat.g.ex.b[1:((k - 1) * (n.eq)), ] <- 0                             #force filling the GexB matrix with 0 before and original values after k

    mat.z <- t(cbind(mat.g.ex, mat.g.ex.b))                             #bind the regressor and breaking regressor matrices together

    l.beta.sigma <- Beta(mat.z = mat.z
                         , mat.y.ex = mat.y.ex
                         , n.eq = n.eq
                         , p
                         , est.mode = est.mode
                         , iter = iter)                                  #compute the BetaSigma object list
    
    # lbetaSigma <- compute_beta_one(Z, Yex, nEq, p, estMode, iter)
    
    mat.beta <- l.beta.sigma$mat.beta                                               #get the vector of betas
    mat.sigma <- l.beta.sigma$mat.sigma                                             #get the covariance matrix of errors
    p.beta <- length(mat.beta)                                                      #get the length of the vector of betas

    #create a selection matrix to get only the betas of interest (breakings)
    if(intercept){
                                                                                    #1 - case where only shift in intercept
      mat.r <- matrix(data = 0, nrow = n.eq, ncol = p.beta)
      mat.r[, (p.beta - n.eq + 1):p.beta] <- diag(n.eq)
    }

    if(!intercept){                                                       
                                                                                    #2 - case where all parameters break
      mat.r <- matrix(data = 0, nrow = p.beta / 2 , ncol = p.beta)
      mat.r[, (p.beta / 2 + 1):p.beta] <- diag(p.beta / 2)
    }

    f.stat[k] <- Fstat(mat.r = mat.r                                     #compute the F-statistic for the current k
                       , mat.beta = mat.beta
                       , mat.z = mat.z
                       , p = p
                       , mat.sigma = mat.sigma)
                                                                         
    mat.ci[k, ] <- ConfidenceInterval(mat.g = mat.g
                                      , mat.s = mat.s
                                      , mat.sigma = mat.sigma
                                      , mat.r = mat.r
                                      , mat.beta = mat.beta
                                      , cv = cv
                                      , p = p)                   #compute the confidence interval for the current k

    mean.shift[k] <- mean(mat.r %*% mat.beta)                                       #get the mean intercept shift
  }

  if(pos.break)                                                           #if posBreak is TRUE, limit to positive break detection
    f.stat[mean.shift < 0] <- 0

  # dev.new()
  plot(f.stat)

  g1 <- PlotStats(my.dates = my.dates
                  , my.vars = my.vars
                  , f.stat = f.stat
                  , mat.ci = mat.ci
                  , mat.y = mat.y)
  
  break.ind <- which.max(f.stat)
  break.date <- my.dates[break.ind]
  break.ci <- mat.ci[break.ind, ]
  rownames(mat.y) <- my.dates
  mat.g.t <- t(mat.g)
  rownames(mat.g.t) <- my.dates
  mean.shift <- mean.shift[break.ind]
  max.f <- max(f.stat, na.rm = T)
  trim.dates <- matrix(data = c(my.dates[start.ind], my.dates[end.ind]), nrow = 1, ncol = 2)
  colnames(trim.dates) <- c("begin trim date", "end trim date")
  
  return(list(f.stat = f.stat
              , max.f = max.f
              , conf.interval = break.ci
              , critical.values = cv
              , break.date = break.date
              , mat.y = data.frame(mat.y)
              , mat.g = data.frame(mat.g.t)
              , break.ind = break.ind
              , mean.shift = mean.shift
              , aic.bic = q.opt
              , g1 = g1
              , trim.dates = trim.dates))
}
