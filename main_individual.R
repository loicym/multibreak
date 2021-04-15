#main entry function to compute BLS (1998) test of structural break for all individual series separately.

MainIndividual <- function(mat.y                           #Y is a matrix or vector which will be lagged by (q) to compute a VAR(q)
                 , mat.x = NULL                             #X is a matrix of (contemporaneous) covariates
                 , trend = FALSE                        #trend is a boolean indicating whether a trend vector should be added to the VAR
                 , intercept = TRUE                     #intercept is a boolean indicating whether the test applies on the mean shift (TRUE) or all parameters (FALSE)
                 , ci = c(0.9, 0.95, 0.99)              #ci is the vector of confidence intervals (in growing order) to compute based on the CDF of a V distr.
                 , est.mode = "OLS"                      #estMode can take values of "OLS", "FGLS", "IGLS"
                 , iter = 3                             #in the case of "IGLS", the number of iteration "iter" can be specified.
                 , aic.bic.mode = "AIC"                   #AicbicMode can be "AIC" or "BIC" depending on the maximum criterion to select
                 , q.max = 6                             #qMax is the number of lags (from 1 to qMax) tested to determine the AIC / BIC maximum.
                 , trim = 0.15                          #trim is the percentage parameter to start and end the sample analysis
                 , pos.break = FALSE                     #if we want the algorithm to only detect positive breaks
                 , elec.dates
){
  
  library(ggplot2)
  elec.dates <- as.Date(as.character(elec.dates), format = "%d.%m.%Y")
  my.vars <- colnames(mat.y)                                                   #get the variable names
  my.dates <- rownames(mat.y)                                                  #get the matching date vector
  
  l.indiv <- list()
  g1 <- ggplot()
  
  mat.breaks <- data.frame(matrix(data = NA, nrow = dim(mat.y)[1], ncol = dim(mat.y)[2]))
  rownames(mat.breaks) <- my.dates
  colnames(mat.breaks) <- my.vars
    
  f.vec <- rep(NA, dim(mat.y)[2])
  for(j in 1:dim(mat.y)[2]){
    
    y.indiv <- data.matrix(mat.y[, j])
    colnames(y.indiv)[1] <- my.vars[j]
    rownames(y.indiv) <- my.dates
    
    l.indiv[[j]] <- Main(mat.y = y.indiv
                         , mat.x = mat.x
                         , trend = trend
                         , intercept = intercept
                         , ci = ci
                         , est.mode = est.mode
                         , iter = iter
                         , aic.bic.mode = aic.bic.mode
                         , q.max = q.max
                         , trim = trim
                         , pos.break = pos.break)
    
    max.f.date <- which.max(l.indiv[[j]]$f.stat)                                     #get the index when the break occurs
    cis <- round(l.indiv[[j]]$conf.interval)                                       #get the confidence intervals around the break
    max.f.val <- l.indiv[[j]]$max.f
    y.indiv <- l.indiv[[j]]$mat.y
    fin.dates <- rownames(y.indiv)
    
    f.vec[j] <- max.f.val

    mat.breaks <- mat.breaks[rownames(mat.breaks) %in% fin.dates, ]               #get the right matrix size
    
    start.date90 <- max.f.date - cis[1]
    end.date90 <- max.f.date + cis[1]
    
    start.date95 <- max.f.date - cis[2]
    end.date95 <- max.f.date + cis[2]
    
    start.date99 <- max.f.date - cis[3]
    end.date99 <- max.f.date + cis[3]
    
    print(start.date99)
    if(start.date99 < 1){
      start.date99 <- 1
    }
    
    if(end.date99 > dim(mat.breaks)[1]){
      end.date99 <- dim(mat.breaks)[1]
    }
    
    mat.breaks[start.date99:end.date99, j] <- max.f.val
  }
  
  return(f.vec)
  
  my.vars <- my.vars
  p <- dim(mat.breaks)[1]
  n <- dim(mat.breaks)[2]

  f.stat.vals <- colMeans(mat.breaks, na.rm = T)
  mat.breaks <- apply(mat.breaks, 2, as.double)
  
  fin.dates <- data.matrix(fin.dates)
  print(dim(fin.dates))
  print(dim(mat.breaks))
  print(fin.dates)
  
  mat.breaks <- data.frame(cbind(fin.dates, mat.breaks))

  colnames(mat.breaks)[1] <- "Date"
  colnames(mat.breaks)[2:(n + 1)] <- my.vars
  
  mat.breaks[, 1] <- as.Date(fin.dates, format = "%d.%m.%Y")
  
  mat.breaks <- melt(mat.breaks, id.var = "Date")                                   #reshape to long format
  names(mat.breaks)[2] = "Variables"

  
  mat.breaks$value <- as.double(mat.breaks$value)
  print(mat.breaks)
  g1 <- ggplot(mat.breaks, aes(x = Date, y = value, group = Variables, colour = Variables))
  g1 <- g1  + geom_line(size = 1.5)
  
  
  g1 <- g1 + theme_bw()
  g1 <- g1 + ggtitle("") + xlab("Date") + ylab("F-statistic")
  # g1 <- g1 + scale_y_continuous(expand = c(0, 0))                                                    #force the y axis to start at zero
  g1 <- g1 + scale_x_date(breaks = scales::pretty_breaks(n = 10))

  d <- data.frame(date = elec.dates, value = f.stat.vals)
  print(d)
  
  # g1 <- g1 + geom_vline(data = d, mapping = aes(xintercept = date), color = "black", size = 1)
  g1 <- g1 + geom_point(data = d, mapping = aes(x = date, y = value), size = 4)
  # g1 <- g1 + geom_point(data = matBreaks, size = 4)
  
  return(g1)
}
