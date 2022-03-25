
#' @importFrom numDeriv jacobian
#' @importFrom memisc getSummary
#' @import stats
#' @export
getSummary.lciv <- function(obj, alpha = 0.05, ...){
  s       <- summary(obj)
  cf      <- s$estimate
  Q       <- obj$Q
  y.var   <- obj$y.var
  end.var <- obj$end.var
  # Obtain shares
  shares   <- make.shares(coef(obj), obj)
  jac      <- numDeriv::jacobian(make.shares, coef(obj), obj = obj)
  se       <- sqrt(diag(jac %*% vcov(obj) %*% t(jac)))
  z        <- shares / se
  p        <- 2 * pnorm(-abs(z))
  shares.T <- cbind(shares, se, z, p)
  rownames(shares.T) <- paste("class", 1:Q, "share", sep = ".")
  # Obtain rhos
  rhos     <- make.rhos(coef(obj), obj)
  jac      <- numDeriv::jacobian(make.rhos, coef(obj), obj = obj)
  se       <- sqrt(diag(jac %*% vcov(obj) %*% t(jac)))
  z        <- rhos / se
  p        <- 2 * pnorm(-abs(z))
  shares.R <- cbind(rhos, se, z, p)
  rownames(shares.R) <- paste("class", 1:Q, "rho", sep = ".")
  # Obtain sigmas
  sigmas   <- make.sigmas(coef(obj), obj)
  jac      <- numDeriv::jacobian(make.sigmas, coef(obj), obj = obj)
  se       <- sqrt(diag(jac %*% vcov(obj) %*% t(jac)))
  z        <- sigmas / se
  p        <- 2 * pnorm(-abs(z))
  shares.S <- cbind(sigmas, se, z, p)
  names.s <- c()
  for (q in 1:Q) names.s <- c(names.s, paste("class", q, paste("eq", 1:2, "sigma", sep = "."), sep = "."))
  rownames(shares.S) <- names.s
  # Wlate
  wlate    <- make.wlate(coef(obj), obj)
  #names(wlate) <- paste("class", 1, "wlate", sep = ".")
  jac      <- numDeriv::jacobian(make.wlate, coef(obj), obj = obj)
  se       <- sqrt(diag(jac %*% vcov(obj) %*% t(jac)))
  z        <- wlate / se
  p        <- 2 * pnorm(-abs(z))
  shares.W <- cbind(wlate, se, z, p)
  rownames(shares.W) <- paste("class", 1, "wate", sep = ".")
  # Make table
  cf       <- rbind(cf, shares.R, shares.S, shares.T, shares.W)
  names.pi <- colnames(model.matrix(obj$formula, data = obj$mf, rhs = 3))
  names.H  <- paste("pi", names.pi, sep = ".")
  cval    <- qnorm(1 - alpha/2)
  cf      <- cbind(cf, cf[, 1] - cval * cf[, 2], cf[, 1] + cval * cf[, 2])

  names.x  <- paste("eq.1",   colnames(model.matrix(obj$formula, data = obj$mf, rhs = 1)), sep = ".")
  names.z  <- paste("eq.2",   colnames(model.matrix(obj$formula, data = obj$mf, rhs = 2)), sep = ".")
  names.pi <- colnames(model.matrix(obj$formula, data = obj$mf, rhs = 3))
  names.H <- paste("pi", names.pi, sep = ".")
  all.vars <- unique(c(names.x, names.z, names.H,
                       paste('eq', 1:2, "lnsigma", sep = '.'),
                       "athrho",
                       "rho",
                       paste('eq', 1:2, "sigma", sep = '.'),
                       "share",
                       "wate"))
  # Table by class
  coef <- array(dim = c(length(all.vars),6, Q),
                dimnames = list(all.vars, c("est", "se", "stat", "p", "lwr", "upr"), paste("class", 1:Q, sep = " ")))
  for (q in 1:Q) {
    temp.vars    <- grep(paste0("class.", q), rownames(cf), value = TRUE)
    temp.vars.no <- gsub(paste0("class.", q, "."), "",  temp.vars)
    coef[rownames(coef) %in% temp.vars.no, , q] <- cf[rownames(cf) %in% temp.vars, ]
  }
  # Statistics
  sumstat <- c(logLik = logLik(obj), deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = nrow(obj$gradientObs),
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts, xlevels = obj$xlevels, call = obj$call)
}
