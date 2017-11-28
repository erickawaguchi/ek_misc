library(survival)
library(aftgee)

syntheticBuckleyJames <- function(timevar, delta, X, n.iter = 4, tol = 1E-8, B = 100) {
  #- Create synthetic data using Koul's method
  fit.km <- survfit(Surv(timevar, (1 - delta)) ~ 1)
  s.est  <- fit.km$surv[fit.km$n.censor == 0]
  s.time <- fit.km$time[fit.km$n.censor == 0]
  G      <- stepfun(s.time, c(1, s.est)) #Censoring distribution 
  t.koul <- (timevar * delta) / G(timevar)
  
  my.ctrl <- aftgee.control(maxiter = 1, reltol = tol)
  #- Initial values using linear regression 
  fit.init <- lm(t.koul ~ X)
  b0       <- fit.init$coefficients
  cov0     <- summary(fit.init)$cov.unscaled

  #- Synthetic Buckley-James
  coef.matrix <- matrix(NA, nrow = n.iter + 1, ncol = length(b0))
  var.matrix  <- matrix(NA, nrow = n.iter + 1, ncol = length(b0))
  rownames(coef.matrix) <- rownames(var.matrix)  <- paste0("Iter ", 0:n.iter)
  colnames(coef.matrix) <- colnames(var.matrix)  <- c("Intercept", paste0("V", 1:ncol(X)))
  
  coef.matrix[1, ] <- b0
  var.matrix[1, ]  <- diag(cov0)
  for(i in 1:n.iter) {
    fit.gee <- aftgee(Surv(timevar, delta) ~ X, binit = b0, B = B, corstr = "indep",
                      control = my.ctrl)
    coef.matrix[i + 1, ] <- fit.gee$coefficients[, 2]
    var.matrix[i + 1, ]  <- diag(fit.gee$var.res)
    b0  <- fit.gee$coefficients[, 2]
    eps <- sum((coef.matrix[i + 1, ] -  coef.matrix[i, ])^2)
  }
  
  out <- list()
  out$final.coef   <- coef.matrix[n.iter + 1, ]
  out$final.var    <- fit.gee$var.res
  out$coefficients <- coef.matrix
  out$var          <- var.matrix
  out$n.iter       <- n.iter
  out$eps          <- eps
  
  return(out)
}

