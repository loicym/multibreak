compute_plot_stats <- function(myDates, Variables, fstat, CI, Y)     #function to compute summary statistics and deliver plots
{
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(stargazer)
  
  myDates <- as.Date(myDates, format = "%d.%m.%Y")
  
  
  maxF <- which.max(fstat)                                         #get the index when the break occurs
  cis <- round(CI[maxF, ])                                         #get the confidence intervals around the break
  
                                                                   #hard code for three different confidence intervals
  
  # maxF <- 60
  
  startDate90 <- myDates[maxF - cis[1]]
  endDate90 <- myDates[maxF + cis[1]]

  startDate95 <- myDates[maxF - cis[2]]
  endDate95 <- myDates[maxF + cis[2]]

  startDate99 <- myDates[maxF - cis[3]]
  endDate99 <- myDates[maxF + cis[3]]
  
  p <- dim(Y)[1]
  n <- dim(Y)[2]
  Y <- Y * 100                                                     #get values in percentage
  
  Y <- apply(Y, 2, as.double)
  Y <- data.frame(cbind(myDates, Y))
  
  colnames(Y)[1] <- "Date"
  colnames(Y)[2:(n + 1)] <- Variables
 
  Y[,1] <- myDates
 
  Y <- melt(Y, id.var = "Date")                                   #reshape to long format
  names(Y)[2] = "Variables"
  
  # dev.new()
  g1 <- ggplot(Y, aes(x = Date, y = value, group = Variables, colour = Variables))
  g1 <- g1  + geom_line()
  
  g1 <- g1 + ggtitle("") + xlab("Date") + ylab("Index investment (%)")
  # g1 <- g1 + theme(axis.title.y = element_text(angle = 0))
  
  #g1 <- g1 + coord_cartesian(ylim = c(0, max(Y$value, na.rm = T)))
  g1 <- g1 + scale_y_continuous(expand = c(0,0))                                                    #force the y axis to start at zero
  
  g1 <- g1 + scale_x_date(breaks = scales::pretty_breaks(n = 10))
  # g1 <- g1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10))
  

  g1 <- g1 + theme_bw()
  # g1 <- g1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   # panel.background = element_blank(), axis.line = element_line(colour = "black"))
                                                                                                    #add shaded area for various ci
  g1 <- g1+annotate("rect", xmin = startDate90, xmax = endDate90, ymin = -Inf, ymax = Inf,
                   alpha = .6)
  g1 <- g1+annotate("rect", xmin = startDate95, xmax = endDate95, ymin = -Inf, ymax = Inf,
                    alpha = .4)
  g1 <- g1+annotate("rect", xmin = startDate99, xmax = endDate99, ymin = -Inf, ymax = Inf,
                    alpha = .2)
  
  d <- data.frame(date = myDates[maxF], event = "index investment break")
  print(d)
  # d=data.frame(date=as.Date(c("2004-02-28", "2006-01-31")), event=c("index investment shock", "arbitrage capital shock"))
  
  g1 <- g1 + geom_vline(data = d, mapping = aes(xintercept = date), color = "black", size = 1)      #add the break line
  
  return(g1)
}