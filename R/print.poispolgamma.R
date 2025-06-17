#' @export
print.poispolgamma <- function(model){
 
  p <- ncol(model$beta)
  if(is.null(p)) p <- 0
  
  # Credibility intervals:
  cat("Credibility Intervals for beta:\n")
  print1 <- apply(model$beta, 2, function(x) quantile(x, c(0.025, 0.975)))
  colnames(print1) <- paste0("beta ", 0:(p-1))
  print(print1)
  # Sigma^2, Tau^2, phi
  cat("Credibility Intervals for sigma2, nugget and phi:\n")
  print2 <- apply(model$theta[,1:3], 2, function(x) quantile(x, c(0.025, 0.975)))
  colnames(print2) <- c("sigma^2", "tau^2", "phi")
  print(print2)
  
  if(!is.null(model$alpha)){
    cat(paste0("Average acceptance rate for phi: "), round(mean(model$alpha), 4))
  }
   
}