#main entry function to compute BLS (1998) test of structural break for all individual series separately.

main_individual <- function(Y                           #Y is a matrix or vector which will be lagged by (q) to compute a VAR(q)
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
                 , elecDates
)
{
  elecDates <- as.Date(as.character(elecDates), format = "%d.%m.%Y")
  myVars <- colnames(Y)                                                   #get the variable names
  myDates <- rownames(Y)                                                  #get the matching date vector
  
  lIndiv <- list()
  g1 <- ggplot()
  
  matBreaks <- data.frame(matrix(data = NA, nrow = dim(Y)[1], ncol = dim(Y)[2]))
  rownames(matBreaks) <- myDates
  colnames(matBreaks) <- myVars
    
  for(j in 1:dim(Y)[2])
  {
    Yindiv <- data.matrix(Y[,j])
    colnames(Yindiv)[1] <- myVars[j]
    rownames(Yindiv) <- myDates
    
    lIndiv[[j]] <- main(Yindiv, X, trend, intercept, ci, estMode, iter, aicbicMode, qMax, trim, posBreak)
    
    maxFDate <- which.max(lIndiv[[j]]$fstat)                                     #get the index when the break occurs
    cis <- round(lIndiv[[j]]$confInterval)                                       #get the confidence intervals around the break
    maxFVal <- lIndiv[[j]]$maxF
    Yindiv <- lIndiv[[j]]$Y
    finDates <- rownames(Yindiv)
    
    # print(rownames(matBreaks) %in% finDates)
    # return()
    matBreaks <- matBreaks[rownames(matBreaks) %in% finDates,]               #get the right matrix size
    
    startDate90 <- maxFDate - cis[1]
    endDate90 <- maxFDate + cis[1]
    
    startDate95 <- maxFDate - cis[2]
    endDate95 <- maxFDate + cis[2]
    
    startDate99 <- maxFDate - cis[3]
    endDate99 <- maxFDate + cis[3]
    
    print(startDate99)
    if(startDate99 < 1)
      startDate99 <- 1
    print(endDate99)
    matBreaks[startDate99:endDate99, j] <- maxFVal
    
    
  }
  
  Variables <- myVars
  p <- dim(matBreaks)[1]
  n <- dim(matBreaks)[2]

  FstatVals <- colMeans(matBreaks, na.rm=T)
  
  matBreaks <- apply(matBreaks, 2, as.double)
  matBreaks <- data.frame(cbind(finDates, matBreaks))
  
  print(matBreaks)
  colnames(matBreaks)[1] <- "Date"
  colnames(matBreaks)[2:(n + 1)] <- Variables
  
  matBreaks[,1] <- as.Date(finDates, format = "%d.%m.%Y")
  
  matBreaks <- melt(matBreaks, id.var = "Date")                                   #reshape to long format
  names(matBreaks)[2] = "Variables"
  

  
  matBreaks$value <- as.double(matBreaks$value)
  print(matBreaks)
  g1 <- ggplot(matBreaks, aes(x = Date, y = value, group = Variables, colour = Variables))
  g1 <- g1  + geom_line(size = 1.5)
  
  
  g1 <- g1 + theme_bw()
  g1 <- g1 + ggtitle("") + xlab("Date") + ylab("F-statistic")
  # g1 <- g1 + scale_y_continuous(expand = c(0, 0))                                                    #force the y axis to start at zero
  g1 <- g1 + scale_x_date(breaks = scales::pretty_breaks(n = 10))
  # dev.new()
  
  # matBreaks$value <- matBreaks$value + 10
  
  d <- data.frame(date = elecDates, value = FstatVals)
  print(d)
  
  # g1 <- g1 + geom_vline(data = d, mapping = aes(xintercept = date), color = "black", size = 1)
  g1 <- g1 + geom_point(data = d, mapping = aes(x = date, y = value), size = 4)
  # g1 <- g1 + geom_point(data = matBreaks, size = 4)
  
  # print(matBreaks)
  
  
  return(g1)
  return(lIndiv)
}
