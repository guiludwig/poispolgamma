#' @export
plot.poispolgamma <- function(model){

  full.posterior.beta <- model$beta
  full.posterior.kappa <- model$kappa
  full.posterior.theta <- model$theta
  
  # pdf("posterior.pdf")
  layout(matrix(1:3, nrow = 3))
  matplot(full.posterior.beta[, -1, drop = FALSE], 
          type = "l", lty = 1, lwd = 1,
          ylab = "posterior", xlab = "iteration", main = "Regression coefficients")
  # legend("bottomleft", lty = 1, col = 1:2, 
  #        legend = c("beta1 (true = 0.3)", "beta2 (true = 0.7)"), 
  #        inset = .01, bg = "white")
  matplot(full.posterior.theta[, 1, drop = FALSE], 
          type = "l", lty = 1, lwd = 1,
          ylab = "posterior", xlab = "iteration", main = "Spatial effect variance")
  # legend("topleft", lty = 1, col = 1:3, 
  #        legend = c("sigma.sq (true = 1)", "nugget (true = 0)"), 
  #        inset = .01, bg = "white")
  matplot(full.posterior.theta[, 3, drop = FALSE], 
          type = "l", lty = 1, lwd = 1,
          ylab = "posterior", xlab = "iteration", col = 3, main = "Spatial effect dependence")
  # legend("topleft", lty = 1, col = 3, 
  #        legend = c("phi (true = 2.5)"), 
  #        inset = .01, bg = "white")
  # dev.off()
  
}