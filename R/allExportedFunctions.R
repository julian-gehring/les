######################################################################
## exported functions for les package
######################################################################

######################################################################
## triangWeight
######################################################################
triangWeight <- function(distance, win)  {

  weight <- (1 - abs(distance)/win)/win

  return(weight)
}


######################################################################
## gaussWeight
######################################################################
gaussWeight <- function(distance, win)  {

  weight <- stats::dnorm(distance, sd=win/2)

  return(weight)
}


######################################################################
## rectangWeight
######################################################################
rectangWeight <- function(distance, win)  {

  weight <- rep(1/(2*win), length(distance))

  return(weight)
}


######################################################################
## epWeight
######################################################################
epWeight <- function(distance, win)  {
  
  weight <- 0.75/win*(1 - (distance/win)^2)

  return(weight)
}


######################################################################
## quartWeight
######################################################################
quartWeight <- function(distance, win)  {
  
  weight <- 15/16/win*(1 - (distance/win)^2)^2

  return(weight)
}


######################################################################
## tricubeWeight
######################################################################
tricubeWeight <- function(distance, win)  {
  
  weight <- 35/32/win*(1 - (distance/win)^2)^3

  return(weight)
}
