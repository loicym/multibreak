## n is number of (time) observations
## p is number of variables
## intensity is the intercept size (absolute value)
## when.break is the percentage of time observations at which the break occurs.
## returns a matrix of brownian with increment of mean zero and variance one for the p variables.

Simul <- function(n = 100
                   , p = 5
                   , intensity = 1
                   , when.break = 0.35){
  
  set.seed(123)                                       #optional
  
  pre.break.n <- round(n * when.break)
  post.break.n <- round(n * (1 - when.break))
  
  y.mat <- rbind(matrix(rnorm(pre.break.n * p) + 0
                        , nrow = pre.break.n, ncol = p)
                        , matrix(rnorm(post.break.n * p) + 0 + intensity
                        , nrow = post.break.n, ncol = p))
  
  start.date <- as.Date("10.03.2021", format = "%d.%m.%Y")                                            #date when I wrote this code
  
  simul.dates <- 0:(n - 1) + start.date
  simul.dates <- format(simul.dates, format = "%d.%m.%Y")

  rownames(y.mat) <- simul.dates
  colnames(y.mat) <- LETTERS[1:p]
  
  print(head(y.mat))
  return(y.mat)
} 