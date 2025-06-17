#' @export
predict.poispolgamma <- function(model, newcoords, newcovariates, 
                                 burnin = 0.2, downsample = 0.01){
  
  coordinates <- model$coordinates
  cov.model <- model$cov.model
  p <- model$ncov
  n <- model$nsite
  NGIBBS <- model$ngibbs
  burnin <- round(nrow(model$theta)*burnin)
  newcovariates <- as.matrix(newcovariates) 
  
  if(burnin > 0){
    bb <- round(seq(burnin, nrow(model$theta), length = round(downsample*(nrow(model$theta) - burnin))))
    
    full.posterior.beta <- model$beta[bb, ]
    full.posterior.epsilon <- model$epsilon[bb, ]
    full.posterior.theta <- model$theta[bb, ]
  } else {
    full.posterior.beta <- model$beta
    full.posterior.epsilon <- model$epsilon
    full.posterior.theta <- model$theta
  }

  epsilon <- matrix(0, ncol = nrow(newcoords), nrow = nrow(full.posterior.epsilon))
  predictive <- epsilon # Same dimensions
  allcoords <- rbind(coordinates, newcoords)
  
  for(j in 1:nrow(predictive)){
    BigSigma <- varcov.spatial(allcoords, 
                               cov.model = cov.model, nugget = 0,
                               cov.pars = c(full.posterior.theta[j,1], 
                                            1/full.posterior.theta[j,3]))$varcov
    # Add nugget only to sampled data
    diag(BigSigma)[1:nrow(coordinates)] <- diag(BigSigma)[1:nrow(coordinates)] + rep(full.posterior.theta[j,2], nrow(coordinates))
    # r <- solve(BigSigma[1:nrow(coordinates), 1:nrow(coordinates)], full.posterior.epsilon[j,]) # instead of y - covariates%*%pm.beta 
    # epsilon[j,] <- BigSigma[(nrow(coordinates)+1):(nrow(allcoords)), 1:nrow(coordinates)] %*% r
    #!# change to ordinary kriging might add stability?
    r <- solve(BigSigma[1:nrow(coordinates), 1:nrow(coordinates)], full.posterior.epsilon[j,] - mean(full.posterior.epsilon[j,])) # instead of y - covariates%*%pm.beta 
    epsilon[j,] <- mean(full.posterior.epsilon[j,]) + BigSigma[(nrow(coordinates)+1):(nrow(allcoords)), 1:nrow(coordinates)] %*% r
    predictive[j,] <- newcovariates%*%full.posterior.beta[j,] + epsilon[j,]
  }
  
  # Do for every spatial effect, produces random sample of y'|y
  
  return(list(predictive = predictive,
              epsilon = epsilon))
    
}