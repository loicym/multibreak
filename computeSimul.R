compute_simul <- function(n, p, intensity = 1, whenbreak = 0.35)
{
  set.seed(123)                                     #optional
  
  prebreakN <- round(n * whenbreak)
  postbreakN <- round(n * (1 - whenbreak))
  
  Y <- rbind(matrix(rnorm(prebreakN * p) + 5, nrow = prebreakN, ncol = p), matrix(rnorm(postbreakN * p) + 5 + intensity, nrow = postbreakN, ncol = p))
  startDate <- as.Date("22.01.2020", format = "%d.%m.%Y")                                            #date when I wrote this code
  
  simulDates <- 0:(n - 1) + startDate
  simulDates <- format(simulDates, format = "%d.%m.%Y")

  rownames(Y) <- simulDates
  colnames(Y) <- LETTERS[1:p]
  
  print(head(Y))
  return(Y)
}