#compute the AIC and BIC criteria for lags from 1 to qMax

AicBic <- function(mat.y
                   , q.max
                   , mat.x
                   , trend
                   , intercept){           

  library(stats)                                                   #load stats package
  aic.bic <- matrix(data <- NA, nrow = 2, ncol = q.max)               #create empty matrix for the AIC / BIC criteria
  
  for(q in 1:q.max){
    
    print(paste0("Testing lags number : ", q))
    
    l.conf.matrix <- ConformableMatrix(mat.y, q, mat.x, trend, intercept) #create a list of conformed objects for the estimation
    mat.y.ex <- l.conf.matrix$mat.y.ex
    mat.g.ex <- l.conf.matrix$mat.g.ex
    
    mod <- lm(mat.y.ex~mat.g.ex)                                            #estimate the model with lm
    aic.bic[1, q] <- AIC(mod)                                       #get AIC
    aic.bic[2, q] <- BIC(mod)                                       #get BIC
  }
  
  rownames(aic.bic) <- c("AIC", "BIC")
  colnames(aic.bic) <- paste0("lags = ", 1:q.max)
  return(aic.bic)
}