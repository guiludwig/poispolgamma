#' @export
coef.poispolgamma <- function(model, burnin = 0.2){
  
  N <- nrow(model$theta)
  
  full.posterior.beta <- model$beta[-seq_len(N*burnin), ]
  full.posterior.theta <- model$theta[-seq_len(N*burnin), ]
  
  colnames(full.posterior.beta) <- paste0("beta ", 0:(ncol(full.posterior.beta)-1))
  colnames(full.posterior.theta) <- c("sigma^2", "tau^2", "phi")
  
  return(list(beta = full.posterior.beta,
              theta = full.posterior.theta))
    
}