#' Fit the multivariate count Gibbs sampler
#' 
#' Multivariate Poisson regression Gibbs sampler using the P\'{o}lya-Gamma method
#' 
#' @param y A n x 1 response vector of integer-valued positive numbers. 
#' (to do: matrix case).
#' @param x A n x p matrix of covariates.
#' @param h An integer that controls the negative binomial 
#' approximation to the Poisson.
#' @param starting Initial values of MCMC samplers, with the exception of the 
#' nugget, which is fixed at the value (initial value for spNNGP package)
#' @param priors Only used by spNNGP package
#' @param tuning Variances for MCMC candidates
#' @param useNNGP Whether to use Finley et al.'s spNNGP package to evaluate
#' the latent Gaussian field. Defaults to FALSE.
#' @param useSS Whether to use Spike and Slab priors on beta for variable
#' selection. Defaults to FALSE.
#' 
#' @export
#' @useDynLib poispolgamma
#' @import spNNGP
#' @importFrom Rcpp evalCpp
#' @examples 
#' # Simulated data: Based on Wang and Wang (2018), same as equation (20) but 
#' using space instead of time
#' n <- 100
#' d <- 1 # Can't handle multivariate data yet
#' set.seed(1)
#' coords <- matrix(runif(n*2), ncol = 2) # 200 spatial sites
#' # B0 <- rbind(1:5, c(.3,.1,.5,.4,.6), -c(.05,.03,.07,.01,.1)) # dimensions in columns, p * d
#' b0 <- c(1, .3, .7) # I increase beta_2 a little
#' X <- cbind(rep(1,n), cos(2*pi*coords[,1]), sin(2*pi*coords[,2])) # n * p
#' # Notice geoR's cov.pars is sigma^2 * rho(dist/phi), whereas spNNGP uses sigma^2 * rho(phi*dist)
#' # We will use sigma^2 * rho(phi*dist) for compatibility with spNNGP
#' epsilon <- geoR::grf(grid = coords, cov.pars = c(1, 1/2.5), 
#'                      nsim = d, cov.model = "exponential", 
#'                      nugget = 0, messages = FALSE)
#' # epsilon$data is n x d
#' lambda <- exp(X%*%b0 + epsilon$data) # unidimensional
#' clambda <- as.numeric(lambda) # same as vec(lambda)
#' y <- rpois(length(clambda), clambda)
#' y <- matrix(y, ncol = d)
#' example1 <- poispolgamma.fit(y, x = X, coords = coords, cov.model = "exponential", 
#'                             h = 1000, ngibbs = 100, outer.it = 200, 
#'                             # initial values of phi cannot be at boundary of uniform prior
#'                             starting = list("sigma.sq" = 1, "tau.sq" = 0.01, "phi" = 2.0), # True: 1, 0, 2.5
#'                             tuning = list("phi" = 0.1),
#'                             priors = list("sigma.sq.ig" = c(1, 1), # shape (alpha) and scale (beta), E(pi) = beta/(alpha-1) = scale/(shape-1)
#'                                           "tau.sq.ig" = c(2, 0.05), # shape and scale
#'                                           "phi.unif" = c(1, 20)),
#'                             verboseParam = TRUE) 
#' plot(example1)
#' \dontrun{
#' # Another example: using spNNGP
#' example2 <- poispolgamma.fit(y, x = X, coords = coords, cov.model = "exponential", 
#'                             h = 1000, ngibbs = 100, outer.it = 200, useNNGP = TRUE,
#'                             # initial values of phi cannot be at boundary of uniform prior
#'                             starting = list("sigma.sq" = 1, "tau.sq" = 0.01, "phi" = 2.0), # True: 1, 0, 2.5
#'                             tuning = list("phi" = 1),
#'                             priors = list("sigma.sq.ig" = c(1, 1), # shape (alpha) and scale (beta), E(pi) = beta/(alpha-1) = scale/(shape-1)
#'                                           "tau.sq.ig" = c(2, 0.05), # shape and scale
#'                                           "phi.unif" = c(1, 20)),
#'                             verboseParam = TRUE) 
#' plot(example2)
#' # Another example: using Spike and Slab priors
#' example3 <- poispolgamma.fit(y, x = X, coords = coords, cov.model = "exponential", 
#'                             h = 1000, ngibbs = 100, outer.it = 200, useSS = TRUE,
#'                             # initial values of phi cannot be at boundary of uniform prior
#'                             starting = list("sigma.sq" = 1, "tau.sq" = 0.01, "phi" = 2.0), # True: 1, 0, 2.5
#'                             tuning = list("phi" = 1),
#'                             priors = list("sigma.sq.ig" = c(1, 1), # shape (alpha) and scale (beta), E(pi) = beta/(alpha-1) = scale/(shape-1)
#'                                           "tau.sq.ig" = c(2, 0.05), # shape and scale
#'                                           "phi.unif" = c(1, 20)),
#'                             verboseParam = TRUE) 
#' plot(example3)
#' }
poispolgamma.fit <- function(y, x, coords, cov.model = "exponential", 
                            h = 1000, ngibbs = 100, outer.it = 200, useNNGP = FALSE, useSS = FALSE,
                            # initial values of phi cannot be at boundary of uniform prior
                            starting = list("var" = 1, "nugget" = 0.01, "phi" = 2.0), 
                            tuning = list("phi" = 0.1), # Metropolis step
                            priors = list("beta.norm.var" = 100,
                                          "sigma.sq.ig" = c(5, 20), # shape (alpha) and scale (beta), E(pi) = beta/(alpha-1) = scale/(shape-1)
                                          "tau.sq.ig" = c(5, 20), # shape and scale
                                          "phi.unif" = c(1, 20)), 
                            startingEpsilon = NULL,
                            hyperpar = c(a0 = 1, b0 = 1, c0 = 1), 
                            SShyperpar = c(a.SS = 1, b.SS = 1, d0 = 5, e0 =  20),
                            tolerance = 1e-4, nomp = 1, nneigh = 15,
                            verboseParam = TRUE){
  
  # (Initialization)
  newY <- as.numeric(y)
  n <- length(newY)
  if(!any(class(coords) %in% "matrix")){
    coords <- as.matrix(coords)
    message("Converted 'coords' to matrix object, please verify if results are correct.")
  }
  
  p <- ncol(x)
  if(is.null(p)) p <- 1
  newCoords <- matrix(0, nrow(coords)*ncol(y), 2) # Repeated coordinates...
  newX <- matrix(0, nrow(coords)*ncol(y), p)
  for(j in 1:ncol(y)){
    newCoords[1:nrow(y) + (j-1)*nrow(y), ] <- coords
    newX[1:nrow(y) + (j-1)*nrow(y), ] <- x
  }
  a0 <- hyperpar[1]
  b0 <- hyperpar[2]
  c0 <- hyperpar[3]
  a.SS <- SShyperpar[1]
  b.SS <- SShyperpar[2]
  d0 <- SShyperpar[3]
  e0 <- SShyperpar[4] 
  
  mu <- log(newY+1) - log(h) #!# PG should be PG(h+Y, |mu-log(h)|) so I subtract log(h) here
  w <- numeric(n)
  savePolGamma <- matrix(0, nrow = outer.it, ncol = n)
  latentAlpha <- rep(0, ngibbs*outer.it)
  latentPosteriorSigma <- rep(0, ngibbs*outer.it)
  
  ret.beta <- matrix(0, ncol = p, nrow = ngibbs*outer.it)
  ret.epsilon <- matrix(0, ncol = n, nrow = ngibbs*outer.it)
  ret.theta <- matrix(0, ncol = 3, nrow = ngibbs*outer.it)
  ret.SSparameters <- matrix(0, ncol = 1 + p, nrow = ngibbs*outer.it) # Changed here
  ret.SSz <- matrix(0, ncol = p, nrow = ngibbs*outer.it)
  
  # step 1:
  for(i in 1:n){
    w[i] <- rpolgamma(1, h + newY[i], mu[i])
  }
  savePolGamma[1,] <- w
  
  cat(paste0("Iteration 1 of ",outer.it,".\n"))
  
  # step 2:
  if(useNNGP){
    
    names(starting) <- c("sigma.sq", "tau.sq", "phi")
    latentModel <- spNNGP(newY ~ newX + 0, coords = newCoords, 
                          starting = starting, # simple kriging model
                          n.neighbors = nneigh, method = "latent", 
                          tuning = tuning, priors = priors, cov.model = cov.model,
                          n.samples = ngibbs, n.omp.threads = nomp, 
                          verbose = verboseParam) 
    
    beta  <- as.matrix(latentModel$p.beta.samples) # ngibbs x p
    epsilon <- t(as.matrix(latentModel$p.w.samples)) # ngibbs x n
    theta <- as.matrix(latentModel$p.theta.samples) # ngibbs x q
    # nugget is tau.sq, that is, theta[2];
    # spatial variance is sigma.sq, that is, theta[1]
    # phi is in theta[3]
    
    # last sampled value of posterior distribution instead
    pl.beta <- t(tail(beta, 1))
    pl.epsilon <- t(tail(epsilon, 1))
    pl.theta <- t(tail(theta, 1))
    
    # step 3:
    b <- 1/(w + 1/pl.theta[2])
    a <- (newY - h)/2 + (newX%*%pl.beta + pl.epsilon)/pl.theta[2]
    a <- as.numeric(a*b)
    mu <- rnorm(length(a), a + log(h), sqrt(b)) # updated mu here
    # mu <- as.numeric(newX%*%pl.beta + pl.epsilon) - log(h)
    
  } else if(useSS) {
    
    latentModel <- MCMCgaussSpikeSlab(newY ~ newX + 0, coords = newCoords, 
                                      starting = starting, # simple kriging model
                                      startingEpsilon = startingEpsilon, 
                                      startingSSz = NULL, 
                                      startingSSparameters = NULL,
                                      n.neighbors = nneigh, method = "latent", 
                                      tuning = tuning, priors = priors, cov.model = cov.model,
                                      n.samples = ngibbs, n.omp.threads = nomp, 
                                      w = w, h = h, a0 = a0, b0 = b0, c0 = c0,
                                      a.SS = a.SS, b.SS = b.SS, 
                                      d0 = d0, e0 = e0) 
    
    beta  <- as.matrix(latentModel$p.beta.samples) # ngibbs x p
    epsilon <- as.matrix(latentModel$p.epsilon.samples) # ngibbs x n
    theta <- as.matrix(latentModel$p.theta.samples) # ngibbs x q
    SSparameters <- as.matrix(latentModel$p.ssparameters.samples)
    SSz <- as.matrix(latentModel$p.z.samples)
    
    # last sampled value of posterior distribution instead
    pl.beta <- t(tail(beta, 1))
    pl.epsilon <- t(tail(epsilon, 1))
    pl.theta <- t(tail(theta, 1))
    pl.SSparameters <- t(tail(SSparameters, 1))
    pl.SSz  <- t(tail(SSz, 1))
    
    mu <- as.numeric(newX%*%pl.beta + pl.epsilon) - log(h) #!!#
    #!# mu <- as.numeric(newX%*%pl.beta + pl.epsilon)
    #!# mu <- exp(as.numeric(newX%*%pl.beta + pl.epsilon))
    latentAlpha[1:ngibbs] <- latentModel$alpha
    latentPosteriorSigma[1:ngibbs] <-latentModel$posteriorSigmaParams
    
  } else {
    
    latentModel <- MCMCgauss(newY ~ newX + 0, coords = newCoords, 
                             starting = starting, # simple kriging model
                             startingEpsilon = startingEpsilon,
                             n.neighbors = nneigh, method = "latent", 
                             tuning = tuning, priors = priors, cov.model = cov.model,
                             n.samples = ngibbs, n.omp.threads = nomp, 
                             w = w, h = h, a0 = a0, b0 = b0, c0 = c0) 
    
    beta  <- as.matrix(latentModel$p.beta.samples) # ngibbs x p
    epsilon <- t(as.matrix(latentModel$p.w.samples)) # ngibbs x n
    theta <- as.matrix(latentModel$p.theta.samples) # ngibbs x q
    
    # last sampled value of posterior distribution instead
    pl.beta <- t(tail(beta, 1))
    pl.epsilon <- t(tail(epsilon, 1))
    pl.theta <- t(tail(theta, 1))
    
    mu <- as.numeric(newX%*%pl.beta + pl.epsilon) - log(h) #!!#
    #!# mu <- as.numeric(newX%*%pl.beta + pl.epsilon)
    #!# mu <- exp(as.numeric(newX%*%pl.beta + pl.epsilon))
    latentAlpha[1:ngibbs] <- latentModel$alpha
    latentPosteriorSigma[1:ngibbs] <-latentModel$posteriorSigmaParams
    
  }
  
  ret.beta[1:ngibbs, ] <- beta
  ret.epsilon[1:ngibbs, ] <- epsilon
  ret.theta[1:ngibbs, ] <- theta
  if(useSS){
    ret.SSparameters[1:ngibbs, ] <- SSparameters
    ret.SSz[1:ngibbs, ] <- SSz
  }
  
  for(it in 2:outer.it){
    
    # step 1:
    for(i in 1:n){
      w[i] <- rpolgamma(1, h + newY[i], mu[i])
    }
    
    new.starting <- list("var" = pl.theta[1], # sigma^2, spatial variance
                         "nugget" = pl.theta[2], # nugget
                         "phi" = pl.theta[3]) # phi
    
    # step 2:
    if(useNNGP){
      
      names(new.starting) <- c("sigma.sq", "tau.sq", "phi")
      latentModel <- spNNGP(mu ~ newX + 0, coords = newCoords, starting = new.starting, # simple kriging model
                            method = "latent", n.neighbors = nneigh,
                            tuning = tuning, priors = priors, cov.model = cov.model,
                            n.samples = ngibbs, n.omp.threads = nomp,
                            verbose = verboseParam) 
      
      beta  <- as.matrix(latentModel$p.beta.samples) # ngibbs x p
      epsilon <- t(as.matrix(latentModel$p.w.samples)) # ngibbs x n
      theta <- as.matrix(latentModel$p.theta.samples) # ngibbs x q

      # last sampled value of posterior distribution instead
      pl.beta <- t(tail(beta, 1)) # tail preserves dimensions as rows, but I want column vector
      pl.epsilon <- t(tail(epsilon, 1))
      pl.theta <- t(tail(theta, 1))
      
      # step 3:
      b <- 1/(w + 1/pl.theta[2])
      a <- (newY - h)/2 + (newX%*%pl.beta + pl.epsilon)/pl.theta[2]
      a <- as.numeric(a*b)
      mu <- rnorm(length(a), a + log(h), sqrt(b))
      # mu <- as.numeric(newX%*%pl.beta + pl.epsilon) - log(h)
      
    } else if(useSS) {

      if(verboseParam) cat(paste0("Iteration ",it," of ",outer.it,".\n"))
      
      latentModel <- MCMCgaussSpikeSlab(newY ~ newX + 0, coords = newCoords, 
                                        starting = new.starting, # simple kriging model
                                        startingEpsilon = pl.epsilon, 
                                        startingSSz = pl.SSz, 
                                        startingSSparameters = pl.SSparameters,
                                        n.neighbors = nneigh, method = "latent", 
                                        tuning = tuning, priors = priors, cov.model = cov.model,
                                        n.samples = ngibbs, n.omp.threads = nomp, 
                                        w = w, h = h, a0 = a0, b0 = b0, c0 = c0,
                                        a.SS = a.SS, b.SS = b.SS, 
                                        d0 = d0, e0 = e0) 
      
      beta  <- as.matrix(latentModel$p.beta.samples) # ngibbs x p
      epsilon <- as.matrix(latentModel$p.epsilon.samples) # ngibbs x n
      theta <- as.matrix(latentModel$p.theta.samples) # ngibbs x q
      SSparameters <- as.matrix(latentModel$p.ssparameters.samples)
      SSz <- as.matrix(latentModel$p.z.samples)
      
      # last sampled value of posterior distribution instead
      pl.beta <- t(tail(beta, 1))
      pl.epsilon <- t(tail(epsilon, 1))
      pl.theta <- t(tail(theta, 1))
      pl.SSparameters <- t(tail(SSparameters, 1))
      pl.SSz  <- t(tail(SSz, 1))
      
      mu <- as.numeric(newX%*%pl.beta + pl.epsilon) - log(h) 
      latentAlpha[1:ngibbs + (it+1)*ngibbs] <- latentModel$alpha
      latentPosteriorSigma[1:ngibbs + (it+1)*ngibbs] <-latentModel$posteriorSigmaParams
      
    } else {
      
      if(verboseParam) cat(paste0("Iteration ",it," of ",outer.it,".\n"))
      
      latentModel <- MCMCgauss(newY ~ newX + 0, coords = newCoords, 
                               starting = new.starting, # simple kriging model
                               startingEpsilon = pl.epsilon,
                               method = "latent", n.neighbors = nneigh,
                               tuning = tuning, priors = priors, cov.model = cov.model,
                               n.samples = ngibbs, n.omp.threads = nomp,
                               w = w, h = h, a0 = a0, b0 = b0, c0 = c0) 
      
      beta  <- as.matrix(latentModel$p.beta.samples) # ngibbs x p
      epsilon <- t(as.matrix(latentModel$p.w.samples)) # ngibbs x n
      theta <- as.matrix(latentModel$p.theta.samples) # ngibbs x q

      # last sampled value of posterior distribution instead
      pl.beta <- t(tail(beta, 1))
      pl.epsilon <- t(tail(epsilon, 1))
      pl.theta <- t(tail(theta, 1))
      
      mu <- as.numeric(newX%*%pl.beta + pl.epsilon) - log(h)
      latentAlpha[1:ngibbs + (it+1)*ngibbs] <- latentModel$alpha
      latentPosteriorSigma[1:ngibbs + (it+1)*ngibbs] <-latentModel$posteriorSigmaParams
      
    }
    
    ret.beta[1:ngibbs + (it-1)*ngibbs, ] <- beta
    ret.epsilon[1:ngibbs + (it-1)*ngibbs, ] <- epsilon
    ret.theta[1:ngibbs + (it-1)*ngibbs, ] <- theta
    if(useSS){
      ret.SSparameters[1:ngibbs + (it-1)*ngibbs, ] <- SSparameters
      ret.SSz[1:ngibbs + (it-1)*ngibbs, ] <- SSz
    }
    
    savePolGamma[it,] <- w
    
  }
  
  ret2 <- list()
  ret2$beta <- ret.beta
  ret2$epsilon <- ret.epsilon
  ret2$theta <- ret.theta
  ret2$w <- savePolGamma
  ret2$ncov <- p
  ret2$nsite <- n
  ret2$ngibbs <- ngibbs
  ret2$outer.it <- outer.it
  ret2$coordinates <- newCoords
  ret2$cov.model <- cov.model
  ret2$X <- x
  ret2$y <- newY
  
  if(!useNNGP){
    ret2$alpha <- latentAlpha 
    ret2$posteriorSigmaParams <- latentPosteriorSigma
  }
  if(useSS){
    ret2$SSz <- ret.SSz
    ret2$SSparameters <- ret.SSparameters
  }
  
  class(ret2) <- "poispolgamma"
  return(ret2)
  
}
