#main entry function to compute BLS (1998) test of structural break.

main <- function(Y                                      #Y is a matrix or vector which will be lagged by (q) to compute a VAR(q)
                 , X = NULL                             #X is a matrix of (contemporaneous) covariates
                 , trend = FALSE                        #trend is a boolean indicating whether a trend vector should be added to the VAR
                 , intercept = TRUE                     #intercept is a boolean indicating whether the test applies on the mean shift (TRUE) or all parameters (FALSE)
                 , ci = c(0.9, 0.95, 0.99)              #ci is the vector of confidence intervals (in growing order) to compute based on the CDF of a V distr.
                 , estMode = "OLS"                      #estMode can take values of "OLS", "FGLS", "IGLS"
                 , iter = 3                             #in the case of "IGLS", the number of iteration "iter" can be specified.
                 , aicbicMode = "AIC"                   #AicbicMode can be "AIC" or "BIC" depending on the maximum criterion to select
                 , qMax = 6                             #qMax is the number of lags (from 1 to qMax) tested to determine the AIC / BIC maximum.
                 , trim = 0.15                          #trim is the percentage parameter to start and end the sample analysis
                 , posBreak = FALSE                     #if we want the algorithm to only detect positive breaks
)
{
  myVars <- colnames(Y)                                                   #get the variable names

  qOpt <- compute_aicbic(Y, qMax, X, trend, intercept)                    #return the AIC and BIC criteria for lags from 1 to 6
  
  
  print(qOpt)                                                             #print the matrix of AIC and BIC for each lags

  q <- as.numeric(which.max(qOpt[aicbicMode,]))                           #choose the lag q according to the max AIC
  print(paste0("lag with the maximum ", aicbicMode, " is: ", q))

  lconf <- matrix_conformation(Y, q, X, trend, intercept)                 #create a list of conform objects for the estimation

  Yex <- lconf$Yex                                                        #get the conformed (expanded) Yex matrix (for the system, in vector form)
  Gex <- lconf$Gex                                                        #get the conformed G matrix of regressors for the system
  p <- lconf$p                                                            #final number of observations
  G <- lconf$G                                                            #original matrix of regressors
  S <- lconf$S                                                            #selection matrix
  Y <- lconf$Y                                                            #matching original dependent variables matrix
  nEq <- lconf$nEq                                                        #original number of equations / dependent variables
  myDates <- lconf$myDates                                                #matching dates

  print(paste0("The number of equations in the system is: ", nEq))

  fstat <- rep(NA, p)                                                     #create a vector of f_statistics for each k tested
  meanShift <- rep(NA, p)                                                 #create a vector with the evaluated size of the intercept difference
  CI <- matrix(data = NA, nrow = p, ncol = length(ci))                    #create a matrix of confidence intervals for each k tested

  cv <- compute_vdistr_cv(ci)                                             #compute critical values for the vector of confidence intervals proposed

  startInd <- round(trim * p)                                             #start index
  endInd <- round(p - trim * p)                                           #end index

  for(k in startInd:endInd)                                               #loop over the k with a trimming date / burn period
  {
    if(k%%10==0)
      print(paste0("The iteration is at the level k = ", k))            #get an idea of where we are in the loop every 10 iterations
    
    
    GexB <- Gex %*% t(S)
    GexB[1:((k - 1) * (nEq)),] <- 0                                       #force filling the GexB matrix with 0 before and original values after k

    Z <- t(cbind(Gex, GexB))                                              #bind the regressor and breaking regressor matrices together

    lbetaSigma <- compute_beta(Z, Yex, nEq, p, estMode, iter)             #compute the BetaSigma object list
    
    # lbetaSigma <- compute_beta_one(Z, Yex, nEq, p, estMode, iter)
    
    Beta <- lbetaSigma$Beta                                               #get the vector of betas
    Sigma <- lbetaSigma$Sigma                                             #get the covariance matrix of errors
    pBeta <- length(Beta)                                                 #get the length of the vector of betas


    #create a selection matrix to get only the betas of interest (breakings)
    if(intercept)                                                         #1 - case where only shift in intercept
    {
      R <- matrix(data = 0, nrow = nEq, ncol = pBeta)
      R[,(pBeta - nEq + 1):pBeta] <- diag(nEq)
    }

    if(!intercept)                                                        #2 - case where all parameters break
    {
      R <- matrix(data = 0, nrow = pBeta / 2 , ncol = pBeta)
      R[,(pBeta / 2 + 1):pBeta] <- diag(pBeta / 2)
    }

    fstat[k] <- compute_fstat(R, Beta, Z, p, Sigma)                      #compute the F-statistic for the current k
    CI[k, ] <- compute_ci(G, S, Sigma, R, Beta, cv, p)                   #compute the confidence interval for the current k

    meanShift[k] <- mean(R%*%Beta)                                       #get the mean intercept shift
  }

  if(posBreak)                                                           #if posBreak is TRUE, limit to positive break detection
    fstat[meanShift < 0] <- 0

  # dev.new()
  plot(fstat)

  g1 <- compute_plot_stats(myDates, myVars, fstat, CI, Y)
  
  breakInd <- which.max(fstat)
  breakDate <- myDates[breakInd]
  breakCi <- CI[breakInd, ]
  rownames(Y) <- myDates
  Gt <- t(G)
  rownames(Gt) <- myDates
  meanShift <- meanShift[breakInd]
  maxF = max(fstat, na.rm=T)
  trimDates = matrix(data=c(myDates[startInd], myDates[endInd]),nrow = 1, ncol = 2)
  colnames(trimDates) <- c("begin trim date", "end trim date")
  
  return(list(fstat = fstat
              , maxF = maxF
              , confInterval = breakCi
              , criticalValues = cv
              , breakDate = breakDate
              , Y = data.frame(Y)
              , G = data.frame(Gt)
              , breakInd = breakInd
              , meanShift = meanShift
              , aicbic = qOpt
              , g1 = g1
              , trimDates = trimDates))
}
