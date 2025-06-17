#' @export
residuals.poispolgamma <- function(model, offset = 1, posteriorQuartile = 0,
                                   burnin = 0.2, downsample = 0.01, 
                                   names = NULL, horizontal = FALSE, 
                                   type = c("posteriorMean", "posteriorSample", "pearson"), ...){
  
  type <- match.arg(type)
  coordinates <- model$coordinates
  cov.model <- model$cov.model
  p <- model$ncov
  n <- model$nsite
  NGIBBS <- model$ngibbs
  burnin <- round(nrow(model$theta)*burnin)
  X <- model$X
  Y <- model$y
  
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
  
  me <- pmax(exp(full.posterior.beta %*% t(X) + full.posterior.epsilon - log(offset)), 0)
  
  if(!is.null(names)){
    colnames(me) <- names
  }
  
  if(posteriorQuartile == 0 | posteriorQuartile == 1){
    mlim <- apply(me, 2, range)
    pc <- ifelse(Y >= mlim[1,] & Y <= mlim[2,], 1, 4)
  } else if(posteriorQuartile > 0 & posteriorQuartile < 0.5) {
    mlim <- apply(me, 2, function(x) quantile(x, c(posteriorQuartile, 1 - posteriorQuartile)))
    pc <- ifelse(Y >= mlim[1,] & Y <= mlim[2,], 1, 4)
  } else if(posteriorQuartile < 1 & posteriorQuartile > 0.5) {
    mlim <- apply(me, 2, function(x) quantile(x, c(1 - posteriorQuartile, posteriorQuartile)))
    pc <- ifelse(Y >= mlim[1,] & Y <= mlim[2,], 1, 4)
  } else {
    stop("posteriorQuantile must be a value in [0,1] except 0.5.")
  }
  
  if(type == "pearson"){
    
    mm <- apply(me, 2, mean)
    ms <- apply(me, 2, sd)
    
    if(horizontal) {
      plot((Y - mm)/ms, 1:length(Y), ...)
    } else {
      plot(1:length(Y), (Y - mm)/ms, ...)
    }
    
    invisible(list(posteriorAverageRate = me, 
                   residualRange = t(mlim),
                   predictedCoverage = mean(pc == 1)))
    
    
  } else if(type == "posteriorMean"){
    
    boxplot(me, ylim = range(c(range(mlim), range(Y))), 
            horizontal = horizontal, ...)
    if(horizontal) {
      points(Y, 1:length(Y), pch = pc, col = "red", cex = 1)
    } else {
      points(Y, pch = pc, col = "red", cex = 1)
    }
    
    invisible(list(posteriorAverageRate = me, 
                   residualRange = t(mlim),
                   predictedCoverage = mean(pc == 1)))
    
  } else {
    
    mtemp <- me
    for(j in 1:ncol(me)){
      for(i in 1:nrow(me)){
        me[i,j] <- rpois(1, mtemp[i,j])
      }
    }
    
    if(!is.null(names)){
      colnames(me) <- names
    }
    
    if(posteriorQuartile == 0 | posteriorQuartile == 1){
      mlim <- apply(me, 2, range)
      pc <- ifelse(Y >= mlim[1,] & Y <= mlim[2,], 1, 4)
    } else if(posteriorQuartile > 0 & posteriorQuartile < 0.5) {
      mlim <- apply(me, 2, function(x) quantile(x, c(posteriorQuartile, 1 - posteriorQuartile)))
      pc <- ifelse(Y >= mlim[1,] & Y <= mlim[2,], 1, 4)
    } else if(posteriorQuartile < 1 & posteriorQuartile > 0.5) {
      mlim <- apply(me, 2, function(x) quantile(x, c(1 - posteriorQuartile, posteriorQuartile)))
      pc <- ifelse(Y >= mlim[1,] & Y <= mlim[2,], 1, 4)
    } else {
      stop("posteriorQuantile must be a value in [0,1] except 0.5.")
    }
    
    boxplot(me, ylim = range(c(range(mlim), range(Y))), 
            horizontal = horizontal, ...)
    if(horizontal) {
      points(Y, 1:length(Y), pch = pc, col = "red", cex = 1)
    } else {
      points(Y, pch = pc, col = "red", cex = 1)
    }
    
    invisible(list(posteriorAverageRate = mtemp, 
                   posteriorSample = me, 
                   residualRange = t(mlim),
                   predictedCoverage = mean(pc == 1)))
    
  }

}