rigamma <- function(n, mu, lambda){
  # Generates random numbers with Inverse-Gamma
  # If X ~ G(mu,lambda) then 1/X ~ IG(mu,lambda) (shape-rate)
  return(1/rgamma(n, mu, lambda))
}