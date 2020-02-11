compute_aicbic <- function(Y, qMax, X, trend, intercept)           #compute the AIC and BIC criteria for lags from 1 to qMax
{
  library(stats)                                                   #load stats package
  AICBIC = matrix(data <- NA, nrow = 2, ncol = qMax)               #create empty matrix for the AIC / BIC criteria
  
  for(q in 1:qMax)
  {
    print(paste0("Testing lags number : ", q))
    
    lConfMatrix <- matrix_conformation(Y, q, X, trend, intercept) #create a list of conformed objects for the estimation
    Yex <- lConfMatrix$Yex
    Gex <- lConfMatrix$Gex
    
    mod <- lm(Yex~Gex)                                            #estimate the model with lm
    AICBIC[1,q] <- AIC(mod)                                       #get AIC
    AICBIC[2,q] <- BIC(mod)                                       #get BIC
  }
  rownames(AICBIC) <- c("AIC", "BIC")
  colnames(AICBIC) <- paste0("lags = ", 1:qMax)
  return(AICBIC)
}