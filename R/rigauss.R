rigauss <- function(n, mu, lambda){
  # Generates random numbers with Inverse-Gaussian distribution
  nu <- rnorm(n)
  y <- nu^2
  x <- mu + (y*mu^2)/(2*lambda) - (mu/(2*lambda))*sqrt(4*mu*lambda*y + (mu*y)^2)
  u <- runif(n)
  return(ifelse(u < mu/(mu+x), x, mu*mu/x))
}