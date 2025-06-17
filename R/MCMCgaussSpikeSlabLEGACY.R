MCMCgaussSpikeSlabLEGACY <- function(formula, data = parent.frame(), coords, starting, 
                                     startingEpsilon = startingEpsilon, 
                                     startingSSz = NULL, startingSSparameters = NULL,
                                     n.neighbors = nneigh, method = "latent", 
                                     tuning = tuning, priors = priors, cov.model = cov.model,
                                     n.samples = ngibbs, n.omp.threads = nomp,
                                     w = w, h = h, a0 = a0, b0 = b0, c0 = c0,
                                     a.SS = a.SS, b.SS = b.SS, alpha1.SS = alpha1.SS, alpha2.SS = alpha2.SS){
  # Arguments match spNNGP arguments, only relevant arguments are 
  # formula, coords, 
  # starting, tuning, priors, cov.model
  # n.samples, verbose
  # phi ~ use metropolis-hastings exponential proposal
  # https://stackoverflow.com/questions/13217322/how-to-reliably-get-dependent-variable-name-from-formula-object
  if (class(formula) != "formula") 
    stop("error: formula is misspecified")
  fm <- terms(formula, data = data, specials = NULL)
  if (missing(data)) 
    data <- sys.frame(sys.parent())
  getY <- get(all.vars(fm)[1], envir = parent.frame())
  getX <- get(all.vars(fm)[2], envir = parent.frame())
  x.names <- attr(fm, "term.labels")
  
  p.beta.samples <- matrix(0, ncol = ncol(getX), nrow = n.samples + 1)
  p.theta.samples <- matrix(0, ncol = 3, nrow = n.samples + 1)
  p.w.samples <- matrix(0, ncol = n.samples + 1, nrow = nrow(getX))
  p.z.samples <- matrix(0, ncol = ncol(getX), nrow = n.samples + 1)
  p.ssparameters.samples <- matrix(0, ncol = 2, nrow = n.samples + 1) # pi, tau^2
  alpha <- numeric(n.samples)
  posteriorSigmaParams <- numeric(n.samples)
  
  # Initial values;  
  if(is.null(startingSSparameters)) {
    #!!!# p.ssparameters.samples[1,] <- c(rbeta(1, a.SS, b.SS), 100) # pi, tau^2 
    p.ssparameters.samples[1,] <- c(rbeta(1, a.SS, b.SS), 1/rgamma(1, alpha1.SS, alpha2.SS)) # pi, tau^2
  } else {
    p.ssparameters.samples[1,] <- startingSSparameters
  }
  if(is.null(startingSSz)){
    p.z.samples[1,] <- c(1, rep(1, ncol(getX)-1))
  } else {
    p.z.samples[1,] <- startingSSz
  }
  # Unless the number of parameters is too large,
  # Computing the inverse is ok
  # Since I'm updating w on blocks, can be done once
  #!!!!# intercept fixed p.ssparameters.samples[1,2]*diag(ncol(getX)) # pi, tau.sq
  #!!!!#  SigmaBeta <- SigmaBeta.inv <- 100 * diag(ncol(getX))
  SigmaBeta <- SigmaBeta.inv <- diag(c(1/100, 1/rep(p.ssparameters.samples[1,2], ncol(getX)-1))) # tau^2 
  A1 <- SigmaBeta.inv + crossprod(getX, diag(w))%*%getX
  A <- solve(A1)
  
  # We can customize prior variance on beta here
  #!# start at previous values if available
  if(is.null(startingEpsilon)){
    p.beta.samples[1,] <- as.numeric(rmvnorm(1, rep(0, ncol(getX)), A)) # Draw from prior
    p.beta.samples[1,] <- p.z.samples[1,]*p.beta.samples[1,] # Same dimensions, will add point mass at 0
  } else {
    a <- -0.5*as.numeric(A%*%crossprod(getX, 2*diag(w)%*%startingEpsilon + h - getY - 2*w*log(h))) #!!# 
    p.beta.samples[1,] <- as.numeric(rmvnorm(1, a, A)) # Draw from prior
    p.beta.samples[1,] <- p.z.samples[1,]*p.beta.samples[1,] # Same dimensions, will add point mass at 0
  }
  
  p.theta.samples[1,] <- c(starting[["var"]], 
                           starting[["nugget"]],
                           starting[["phi"]]) # initial values
  SB <- varcov.spatial(coords, cov.model = cov.model,
                       nugget = p.theta.samples[1,2],
                       cov.pars = c(p.theta.samples[1,1], 
                                    1/p.theta.samples[1,3]),
                       inv = TRUE)$inverse
  B <- solve(SB + diag(w))
  b <- -0.5*B%*%(2*diag(w)%*%getX%*%p.beta.samples[1,] + h - getY - 2*w*log(h)) #!!#
  p.w.samples[,1] <- as.numeric(rmvnorm(1, b, B)) 
  
  for(i in 2:(n.samples+1)){
    
    # Draw Spike and Slab parameters
    # p.ssparameters.samples[i,] <- c(rbeta(1, a.SS + sum(p.z.samples[i-1,]), b.SS + sum(1-p.z.samples[i-1,])), 
    #                                 100) # pi, tau^2
    # No intercept!!
    p.ssparameters.samples[i,] <- c(rbeta(1, a.SS + sum(p.z.samples[i-1,-1]), b.SS + sum(1-p.z.samples[i-1,-1])),
                                    1/rgamma(1, alpha1.SS + 0.5*sum(p.z.samples[i-1,-1]),
                                             alpha2.SS + 0.5*sum(p.z.samples[i-1,-1]*p.beta.samples[i-1,-1]^2))) # pi, tau^2
    # Sample Spike and Slab Z: ############# REVISE HERE
    r.ij <- 0*getX # Same number of rows and columns
    probabilities <- numeric(ncol(getX))
    for(ell in 2:ncol(getX)){ 
      r.ij[,ell] <- as.numeric(getX[,-ell]%*%p.beta.samples[i-1,-ell]) - log(h) 
      c.j.denom <- 1/p.ssparameters.samples[i,2] + sum(as.numeric(w) * getX[,ell]^2) # pi, tau
      c.j.num <- sum(as.numeric(w) * r.ij[,ell] * getX[,ell] - 0.5*(h-getY)*getX[,ell])
      C.j <- sqrt(2*pi/c.j.denom)*exp(0.5*(c.j.num^2)/c.j.denom)
      # Maybe save these probabilities?
      probabilities[ell] <- p.ssparameters.samples[i,1]/(p.ssparameters.samples[i,1] + (1 - p.ssparameters.samples[i,1])/C.j)
    }
    probabilities[1] <- 1 # Intercept never gets penalized

    p.z.samples[i,] <- rbinom(ncol(getX), 1, prob = probabilities)

    # Draw beta
    # Remember p.w.samples has spatial effects in columns
    # SigmaBeta.inv <- (1/100) * diag(ncol(getX))
    #!!!!# intercept fixed, not penalizing
    SigmaBeta.inv <- diag(c(1/100, rep((1/p.ssparameters.samples[i,2]), ncol(getX)-1))) # tau^2 
    A1 <- SigmaBeta.inv + crossprod(getX, diag(w))%*%getX
    A <- solve(A1)
    a <- -0.5*as.numeric(A%*%crossprod(getX, 2*diag(w)%*%p.w.samples[,i-1] + h - getY - 2*w*log(h))) #!!#
    p.beta.samples[i,] <- as.numeric(rmvnorm(1, a, A)) 
    p.beta.samples[i,] <- p.z.samples[i,]*p.beta.samples[i,]
    # Draw epsilon (save as "w" to keep transparent with spNNGP)
    SB <- varcov.spatial(coords, cov.model = cov.model,
                         nugget = p.theta.samples[i-1, 2],
                         cov.pars = c(p.theta.samples[i-1, 1], 
                                      1/p.theta.samples[i-1, 3]),
                         inv = TRUE)$inverse
    B <- solve(SB + diag(w))
    b <- -0.5*B%*%(2*diag(w)%*%getX%*%t(p.beta.samples[i,,drop=FALSE]) + h - getY - 2*w*log(h)) #!!#
    p.w.samples[,i] <- as.numeric(rmvnorm(1, b, B)) 
    
    # Updating after spatial effects
    # Sampling phi' from a log-Gaussian distribution
    # with mu = log(phi) and sd = 
    phi.prime <- rlnorm(1, meanlog = log(p.theta.samples[i-1, 3]), sdlog = tuning[["phi"]]) 
    SB.prime <- varcov.spatial(coords, cov.model = cov.model,
                               nugget = p.theta.samples[i-1, 2],
                               cov.pars = c(p.theta.samples[i-1, 1], 
                                            1/phi.prime),
                               inv = TRUE)$inverse
    # Assuming a prior for phi equal to phi ~ c e^(-c phi)
    # Such that E(phi) = initial
    # SB holds P^{-1}
    # Need this copy to stable-eval ratio of determinants...
    SB2 <- varcov.spatial(coords, cov.model = cov.model,
                          nugget = p.theta.samples[i-1, 2],
                          cov.pars = c(p.theta.samples[i-1, 1], 
                                       1/p.theta.samples[i-1, 3]))$varcov
    # instead of sqrt(det(SB.prime)/det(SB)), note SB
    rpost <- sqrt(det(SB2%*%SB.prime))*exp(-0.5*crossprod(p.w.samples[,i], 
                                                          crossprod(SB.prime-SB, 
                                                                    p.w.samples[,i]))
                                           -c0*(phi.prime-p.theta.samples[i-1, 3])) # *starting[["phi"]]
    # My proposals are Q(phi'|phi[t-1]) ~ logN(mu = log(phi[t-1]), sigma = tuning[["phi"]])
    # Ratio should be Q(phi[t-1]|phi')/Q(phi'|phi[t-1])
    # rproposals <- dlnorm(p.theta.samples[i-1, 3],
    #                      meanlog = log(phi.prime),
    #                      sdlog = tuning[["phi"]])/dlnorm(phi.prime, 
    #                                                      meanlog = log(p.theta.samples[i-1, 3]), 
    #                                                      sdlog = tuning[["phi"]])
    rproposals <- phi.prime/p.theta.samples[i-1, 3]
    alpha.prime <- min(1, rpost*rproposals)
    if(is.na(alpha.prime)){
      alpha.prime <- 0
    } # Won't happen if n = 100, but may happen if n = 400
    if(runif(1) > alpha.prime){ # if U < alpha.prime, accept proposal
      # So I am just considering the reject case here
      phi.prime <- p.theta.samples[i-1,3] # previous value
    }
    # To remove the spatial variance from the inverse, so it is a spatial 
    # correlation matrix now.
    PB <- varcov.spatial(coords, cov.model = cov.model,
                         nugget = p.theta.samples[i-1, 2],
                         cov.pars = c(1, # Change done here
                                      1/phi.prime),
                         inv = TRUE)$inverse
    # save posterior parameter values
    # Here i-1 because I did not evaluate this for the starting step,
    # which uses values from previous MCMCgauss call. I organize indexes
    # at the return step
    posteriorSigmaParams[i-1] <- 0.5*crossprod(p.w.samples[,i], 
                                               crossprod(PB, 
                                                         p.w.samples[,i]))
    sigma.prime <- rigamma(1,
                           a0 + length(w)/2, # a0 + m/2
                           b0 + posteriorSigmaParams[i-1])
    p.theta.samples[i,] <- c(sigma.prime,
                             starting[["nugget"]],
                             phi.prime) 
    alpha[i-1] <- alpha.prime
    
  }
  
  # For some reason spNNGP returns n x ngibbs 
  # in p.w.samples, we need to transpose in order to match
  # I am not including initial values again since they were part of previous simulation
  ret <- list(p.beta.samples = p.beta.samples[2:(n.samples+1),],
              p.w.samples = p.w.samples[,2:(n.samples+1)],
              p.theta.samples = p.theta.samples[2:(n.samples+1),],
              p.z.samples = p.z.samples[2:(n.samples+1),],
              p.ssparameters.samples = p.ssparameters.samples[2:(n.samples+1),],
              X = getX,
              alpha = alpha,
              posteriorSigmaParams = posteriorSigmaParams)
  
  return(ret)
  
}
