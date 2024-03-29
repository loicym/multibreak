#computes the critical values for a vector of confidence intervals proposed, ci

Vdistr <- function(ci){          

  n <- length(ci)                         #get the number of ci elements
  target <- 1 - (1 - ci) / 2              #redefine target for a two tail CI
  
  print(paste0("The vdistr targets are: ", target))
  x <- seq(-200, 200, 0.01)               #define the support sequence "x" for the CDF of V
  
  mat.v <- (3 / 2) * exp(abs(x)) * pnorm( (-3 / 2) * abs(x)^0.5 )  - (1 / 2) * pnorm( (-1 / 2) * abs(x)^0.5 ) #compute V
  
  cum.v <- cumsum(mat.v) / sum(mat.v)                #scale the CDF of V to reach one
  
  
  # dev.new()                                               
  # plot(x, cumsum(gamma)/sum(gamma), t = 'l')                    #optionally plot V (nice shape!)
  
  cv <- rep(NA, n)
  k <- 1
  
  print(target)
  for(i in 2:length(x)){
    
    if(cum.v[i - 1] < target[k] && cum.v[i] >= target[k]){
      
      cv[k] <- x[i]
      k <- k + 1
      if(k > n)
        break
    }  
  }
  return(cv)
}