#### IV-het with discrete heterogeneity

ivhetLc <- function(formula, data, subset, na.action,
                     method     = "bfgs", 
                     model      =  c("lciv", "lc"),
                     Q          = 2, 
                     gradient   = TRUE,
                     print.init = FALSE, 
                     messages   = TRUE,
                     start      =  NULL,
                     pert.init  = 0.001, 
                     ...){
  # Required packages
  require("Formula", quietly = TRUE)
  require("maxLik",  quietly = TRUE)

  
  callT     <- match.call(expand.dots = TRUE)
  callF     <- match.call(expand.dots = FALSE)
  nframe    <- length(sys.calls())
  model     <- match.arg(model)
  #gradient  <- match.arg(gradient)

  # ============================
  # 2. Initial checks 
  # ============================
  
  if (Q < 2) stop("Classes cannot be lower than 2")
  
  # ============================
  # 3. Model frame 
  # ============================
  mf         <- match.call(expand.dots = FALSE)
  m          <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf         <- mf[c(1, m)]
  f          <- Formula(formula)
  mf[[1]]    <- as.name("model.frame")
  mf$formula <- f
  mf         <- eval(mf, parent.frame())
  if (length(f)[2] < 2L) stop("You need to specify the instruments in the second part of the formula")

  
  # ============================
  # 4. Make variables 
  # ============================
  
  # Extract variables for each part
  if (model == "lciv") {
    y1    <- model.response(mf)                                 # Continuous dependent variable
    y.var <- f[[2]]                                             # Name of continuous dependent variable
    X     <- model.matrix(f, data = mf, rhs = 1)                # Variables in the first equation
    Z     <- model.matrix(f, data = mf, rhs = 2)                # Variables in the second equation
    y2    <- X[, !(colnames(X) %in% colnames(Z)), drop = FALSE] # Endogenous continuous variable
    end.var     <- colnames(y2)                                 # Name of endogenous variable
    instruments <- colnames(Z)                                  # Name of instruments 
    # Checks
    if (ncol(Z) < ncol(X)) stop("Model underidentified: check your variables in formula")
    if (messages && (ncol(Z) == ncol(X))) cat("\nEstimating a just identified model....\n")
    if (messages && (ncol(Z) > ncol(X)))  cat("\nEstimating an overidentified model....\n")
    if (length(end.var) > 1L) stop(" ivhetLc only works with one endogenous variable") 
  } else {
    y1 <- model.response(mf)
    X  <- model.matrix(f, mf)
  }
  
  # ============================
  # 4. Create H matrix 
  # ============================
  H         <- model.matrix(f, data = mf, rhs = 3)
  names.h   <- colnames(H) 
  indata    <- mf
  cldata    <- indata[rep(seq_len(nrow(indata)), each = Q), ] #expand H
  class     <- factor(rep(1:Q, nrow(X)))
  N         <- length(y1)  
  cldata    <- cbind(cldata, class, row.names = NULL)
  class.var <- formula(f, rhs = 3, lhs = 0) 
  has.int   <- attr(terms(class.var), "intercept") == 1L
  if (has.int) intercept.char <- "factor(class)" else intercept.char <- NULL 
  if (!has.int) {
    class.var      <- update(class.var, ~. + 1)
    class.var      <- update(class.var, ~.)
    class.var.char <- as.character(class.var)[2]
    has.xclass     <- as.character(class.var)[2]
    class.var.var  <- colnames(model.matrix(update(class.var, ~. + 1), cldata))[-1]
    class.var.char <- paste("(", class.var.var, "):class", sep = "")
    if (class.var.char == "1") class.var.char <- class.var.var <- NULL
  } else {
    has.xclass <- as.character(class.var)[2]
    if (has.xclass == "1") {
      class.var.char <- NULL
    } else {
      class.var.var  <- colnames(model.matrix(update(class.var, ~. + 1), cldata))[-1]
      class.var.char <- paste("(", class.var.var, "):class", sep = "")
    } 
  }
  resp.name <- as.character(attr(f, "lhs"))
  form.char <- paste(c(intercept.char, class.var.char), collapse = "+")
  nformula  <- as.formula(paste(resp.name, " ~ ", form.char))
  H         <- model.matrix(nformula, cldata)[, -1, drop = F]
  
  lev1 <- levels(class)[1]
  lev1 <- paste("class", lev1, sep = "")
  if (has.xclass != "1") {
    toremove     <- unlist(lapply(as.list(class.var.var), function(x) paste(lev1, x, sep = ":")))
    revtoremove  <- unlist(lapply(as.list(class.var.var), function(x) paste(x, lev1, sep = ":")))
    toremove     <- colnames(H) %in% c(toremove, revtoremove)
    H            <- H[, !toremove, drop = FALSE]
  } 
  #namesH <- colnames(H)
  #for (i in 1:length(namesH)) namesH[i] <- sub('factor', '', namesH[i])
  #colnames(H) <- namesH
  namesH <- c()
  for (i in 1:length(names.h)) namesH <- c(namesH, paste("class", 2:Q, "pi", names.h[i], sep = "."))
  
  
  Hl <- vector(length = Q, mode = "list")
  names(Hl) <- levels(class)
  for (i in levels(class)) {
    Hl[[i]] <- H[class == i, , drop = FALSE]
  }
  
  # ============================
  # 4. Initial values
  # ============================
  if (is.null(start)){
    s_class <- rep(0, ncol(H))
    names(s_class) <- namesH
    if (model == "lc") {
      # Estimate standard OLS and +/- something
      if(messages) cat("Estimating a OLS model for initial values", fill = TRUE)
      fl         <- formula(f, rhs = 1, lhs = 1)
      ols        <- lm(fl, data = mf)
      beta       <- coef(ols)
      init.shift <- seq(-0.05, 0.05, length.out = Q)
      lc.beta    <- c()
      lc.nbeta   <- c()
      for (i in 1:Q) {
        lc.beta  <- c(lc.beta,  beta + init.shift[i])
        lc.nbeta <- c(lc.nbeta , paste('class', i, names(beta) , sep = '.'))
      }
      names(lc.beta)  <- lc.nbeta
      start <- c(lc.beta, s_class)
    } 
     else {
      # Estimate both equations separately 
      if(messages) cat("Estimating two OLS models for initial values", fill = TRUE)
      ols1         <- lm(y1 ~ X - 1)
      ols2         <- lm(y2 ~ Z - 1)
      beta         <- coef(ols1)
      delta        <- coef(ols2)
      lnsigma_e    <- log(sigma(ols1))
      lnsigma_v    <- log(sigma(ols2))
      rho          <- cor(resid(ols1), resid(ols2))
      athrho       <- 0.5 * log((1 + rho)/(1 - rho))
      names(beta)  <- paste("eq.1", colnames(X), sep = ".")
      names(delta) <- paste("eq.2", colnames(Z), sep = ".")
      names(lnsigma_e) <- paste("eq.1", "lnsigma", sep = ".")
      names(lnsigma_v) <- paste("eq.2", "lnsigma", sep = ".")
      init.shift   <- seq(-pert.init, pert.init, length.out = Q)
      lc.beta      <- lc.delta  <- lc.lnsigma_e  <- lc.lnsigma_v  <- lc.athrho  <- c()
      lc.nbeta     <- lc.ndelta  <- lc.nlnsigma_e  <- lc.nlnsigma_v  <- lc.nathrho  <- c()
      for (i in 1:Q) {
        lc.beta       <- c(lc.beta,      beta      + init.shift[i])
        lc.delta      <- c(lc.delta,     delta     + init.shift[i])
        lc.lnsigma_e  <- c(lc.lnsigma_e, lnsigma_e + init.shift[i])
        lc.lnsigma_v  <- c(lc.lnsigma_v, lnsigma_v + init.shift[i])
        lc.athrho     <- c(lc.athrho,    athrho    + init.shift[i])
        lc.nbeta      <- c(lc.nbeta ,     paste('class', i, names(beta)      , sep = '.'))
        lc.ndelta     <- c(lc.ndelta,     paste('class', i, names(delta)     , sep = '.'))
        lc.nlnsigma_e <- c(lc.nlnsigma_e, paste('class', i, names(lnsigma_e) , sep = '.'))
        lc.nlnsigma_v <- c(lc.nlnsigma_v, paste('class', i, names(lnsigma_v) , sep = '.'))
        lc.nathrho    <- c(lc.nathrho,    paste('class', i, "athrho"         , sep = '.'))
      }
      names(lc.beta)      <- lc.nbeta
      names(lc.delta)     <- lc.ndelta
      names(lc.lnsigma_e) <- lc.nlnsigma_e
      names(lc.lnsigma_v) <- lc.nlnsigma_v
      names(lc.athrho)    <- lc.nathrho
      start <- c(lc.beta, lc.delta, s_class, lc.lnsigma_e, lc.lnsigma_v, lc.athrho)
     }
  } 
  else {
      start <- start
} 

 if (print.init) {
    cat("Initial values:", "\n")
    print(start)
  }

  
  # ============================
  # 5. Estimate the model
  # ============================
  opt <- callT
  opt$start <- start
  #opt$gradient <- as.name('gradient')
  m <- match(c('method', 'print.level', 'iterlim',
               'start','tol', 'ftol', 'steptol', 'fixed', 'constraints', 
               'control', 'gradient'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]] <- as.name("maxLik")
  if(messages) cat("Estimating an IVLC-Linear model", "\n")
  opt$logLik <- as.name("ml_lcivd2")
  opt[c("y1", "y2", "X", "Z", "W", "Q")] <- list(as.name("y1"), 
                                                     as.name("y2"), 
                                                     as.name("X"), 
                                                     as.name("Z"), 
                                                     as.name("Hl"), 
                                                     as.name("Q")) 
  #print(opt)
  out <- eval(opt, sys.frame(which = nframe))
  
  ##############################
  ## ave results
  ##############################
  out$y1          <- y1
  out$y2          <- y2
  out$X           <- X
  out$Z           <- Z
  out$Hl          <- Hl
  out$H           <- H
  out$end.var     <- end.var
  out$instruments <- instruments
  out$formula     <- f
  out$mf          <- mf
  out$call        <- callT
  out$model       <- model
  out$Q           <- Q
  out$y.var       <- y.var
  class(out) <- c("lciv", "maxLik", class(out))
  return(out)
}  



ml_lcivd2 <- function(theta, y1, y2, X, Z, W, Q, gradient = TRUE){
  #param: bq, dq, lq, lnsigma_eq, lnsigma_vq, athrho_q
  K <- ncol(X)
  P <- ncol(Z)
  L <- ncol(W[[1]])
  N <- nrow(X)
  
  beta      <- matrix(theta[1L:(K * Q)], nrow = K, ncol = Q)                        # Matrix of K*Q
  delta     <- matrix(theta[((K * Q) + 1):((P * Q) + (K * Q))], nrow = P, ncol = Q) # Matrix of P*Q
  lambda    <- theta[((P * Q) + (K * Q) + 1):(L + (P * Q) + (K * Q))]               # vector of (Q - 1) * L
  lnsigma_e <- repRows(tail(theta, n = 3 * Q)[1:Q], n = nrow(X))                    # Matrix of n * Q
  lnsigma_v <- repRows(tail(theta, n = 2 * Q)[1:Q], n = nrow(X))                    # Matrix of n * Q
  athrho    <- repRows(tail(theta, n = Q), n = nrow(X))                             # Matrix of n * Q
  sigma_e   <- exp(lnsigma_e)
  sigma_v   <- exp(lnsigma_v)
  rho       <- tanh(athrho)
  # To force the optimizer to find a new parameter
  if (any(rho == 1 || rho == -1)) return(NA)
  
  rownames(beta)  <- colnames(X); colnames(beta)  <- paste("class", 1:Q, sep = ":")
  rownames(delta) <- colnames(Z); colnames(delta) <- paste("class", 1:Q, sep = ":")
  
  # Create components for log-likelihood function 
  f      <- dnorm
  index1 <- tcrossprod(X, t(beta))        #N * Q
  index2 <- tcrossprod(Z, t(delta))       #N * Q
  y1     <- repCols(y1, n = Q)            #N * Q
  y2     <- repCols(y2, n = Q)            #N * Q
  a      <- (y1 - index1)
  b      <- (y2 - index2)
  r      <- (sigma_e / sigma_v)
  aiq    <- (a - r * rho * b) / sqrt((1 - rho^2) * sigma_e^2)
  biq    <- b / sigma_v
  #P1     <-  (1 / sqrt((1 - rho^2) * sigma_e^2)) * f(aiq)
  P1     <- dnorm(y1,  mean = index1 + r * rho * b, sd = sqrt((1 - rho^2) * sigma_e^2))
  #P1     <- pmax(P1, .Machine$double.eps)
  #P1[is.nan(P1)] <- 0
  #P2   <- (1 / sigma_v) * f(biq)
  P2     <- dnorm(y2, mean = index2, sd = sigma_v)
  #P1[is.nan(P2)] <- 0
  
  # Make weights from a multinomial logit model
  ew <- lapply(W, function(x) exp(crossprod(t(x), lambda)))
  sew <- suml(ew)
  Wiq <- lapply(ew, function(x){ v <- x / sew;
  v[is.na(v)] <- 0;
  as.vector(v)})
  Wiq <- Reduce(cbind, Wiq) 
  
  Piq <- P1 * P2
  #Piq <- pmax(P1 * P2 , .Machine$double.eps)
  Pi <- apply(Wiq * Piq, 1, sum)
  #Pi <- pmax(apply(Wiq * Piq, 1, sum), .Machine$double.eps)
  Pi <- pmax(Pi, .Machine$double.eps)
  #Pi[is.nan(Pi)] <- Inf
  LL <- log(Pi)

  
  ##Gradient
  if (gradient) {
    frat     <- function(x) (-x * dnorm(x)) / pmax(dnorm(x), .Machine$double.eps)
    dp1  <-  -1 * frat(aiq) * (1 / sqrt((1 - rho^2) * sigma_e^2))
    dp2  <-  (frat(aiq) * (1 / sqrt((1 - rho^2) * sigma_e^2)) * r * rho - frat(biq) * (1 / sigma_v))
    dpse <-  -1 - frat(aiq) * (a / sqrt((1 - rho^2) * sigma_e^2))
    dpsv <-   frat(aiq) * (rho * sigma_e * b / (sqrt((1 - rho^2) * sigma_e^2) * sigma_v)) - frat(biq) * (b / sigma_v) - 1
    dpdt <-   rho + frat(aiq) * ((a * rho - b*r) / sqrt(sech(athrho)^2 * sigma_e^2))
      
    Qiq  <-  (Wiq * Piq / Pi)
    eta  <-  Qiq * dp1  
    etar <- eta[, rep(1:Q, each = K)]
    Xg <- X[, rep(1:K, Q)] # N * (K * Q)
    grad.beta <- Xg * etar
    
    eta <- Qiq * dp2
    etar <- eta[, rep(1:Q, each = P)]
    Zg <- Z[, rep(1:P, Q)] # N * (P * Q)
    grad.delta <- Zg * etar
    
    grad.lnsigma_e <- Qiq * dpse
    grad.lnsigma_v <- Qiq * dpsv
    
    grad.athrho <- Qiq * dpdt
    
    Wg <- vector(mode = "list", length = Q)
    IQ <- diag(Q)
    for (q in 1:Q) Wg[[q]] <- rowSums(Qiq * (repRows(IQ[q, ], N) - repCols(Wiq[, q], Q)))
    grad.lambda <- suml(mapply("*", W, Wg, SIMPLIFY = FALSE))
    grad <- cbind(grad.beta, grad.delta, grad.lambda, grad.lnsigma_e, grad.lnsigma_v, grad.athrho)
    colnames(grad) <- names(theta)
    attr(LL,'gradient') <- grad
  }
  LL
}

sech <- function(x) 1 / cosh(x)

## suml function from mlogit (Croissant, 2013)
suml <- function(x){
  n <- length(x)
  if (!is.null(dim(x[[1]]))) {
    d <- dim(x[[1]])
    s <- matrix(0,d[1], d[2])
    for (i in 1:n) {
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  else{
    s <- rep(0,length(x[[n]]))
    for (i in 1:n) {
      x[[i]][is.na(x[[i]])] <- 0
      s <- s + x[[i]]
    }
  }
  s
}


repRows <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

repCols <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

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
  rownames(shares.W) <- paste("class", 1, "wlate", sep = ".")
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
                       "wlate"))
  # Table by class
  coef <- array(dim = c(length(all.vars),6, Q), 
                dimnames = list(all.vars, c("est", "se", "stat", "p", "lwr", "upr"), paste("class", 1:Q, sep = " ")))
  for (q in 1:Q) {
    temp.vars    <- grep(paste0("class.", q), rownames(cf), value = TRUE)
    temp.vars.no <- gsub(paste0("class.", q, "."), "",  temp.vars)
    coef[rownames(coef) %in% temp.vars.no, , q] <- cf[rownames(cf) %in% temp.vars, ]
  }
  # Statistics
  sumstat <- c(logLik = logLik(obj), deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = nObs(obj), 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts, xlevels = obj$xlevels, call = obj$call)
}


make.shares <- function(theta, obj, ...){
  theta    <- theta
  names.pi <- colnames(model.matrix(obj$formula, data = obj$mf, rhs = 3))
  names.H  <- paste("pi", names.pi, sep = ".")
  H.bar    <- colMeans(model.matrix(obj$formula, data = obj$mf, rhs = 3))
  Q        <- obj$Q
  expsq    <- c(exp(0))
  for (q in 2:Q){
    hd <- c(0)
    for (l in 1:length(H.bar)){
        name.delta <- paste("class", q, names.H[l], sep = ".")
        temp       <- H.bar[l] * theta[name.delta]
        hd         <- hd + temp
    }
    expsq <- c(expsq, exp(hd))
  }
  shares <- expsq / sum(expsq)
  return(shares)
}  

make.rhos<- function(theta, obj, ...){
  Q <- obj$Q
  rhos <- c()
  for (q in 1:Q){
    name.rho <- paste("class", q, "athrho", sep = ".")
    rhos <- c(rhos, tanh(theta[name.rho]))
  }
  return(rhos)
} 

make.sigmas <- function(theta, obj, ...){
  Q      <- obj$Q
  sigmas <- c()
  for (q in 1:Q){
    name.sigma1 <- paste("class", q, "eq.1.lnsigma", sep = ".")
    name.sigma2 <- paste("class", q, "eq.2.lnsigma", sep = ".")
    sigmas      <- c(sigmas, exp(theta[name.sigma1]), exp(theta[name.sigma2]))
  }
  return(sigmas)
}  


make.wlate <- function(theta, obj, ...){
  shares <- make.shares(theta, obj)
  end    <- obj$end.var
  beta.names  <- grep(end, names(theta), value = TRUE)  
  wlate <- sum(shares * theta[beta.names])
  return(wlate)
}  


fit.measures <- function(obj, ...){
    AIC <- AIC(obj)
    BIC <- BIC(obj)
    
}
  

# gen_init <- function(y1, y2, X, Z, H, namesH, Hl, Q, R = 50){
#   s_class        <- rep(0, ncol(H))
#   names(s_class) <- namesH
#   ols1           <- lm(y1 ~ X - 1)
#   ols2           <- lm(y2 ~ Z - 1)
#   beta           <- coef(ols1)
#   delta          <- coef(ols2)
#   lnsigma_e      <- log(sigma(ols1))
#   lnsigma_v      <- log(sigma(ols2))
#   rho            <- cor(resid(ols1), resid(ols2))
#   athrho         <- 0.5 * log((1 + rho)/(1 - rho))
#   
#   #pert   <- runif(R, min = 0, max = 5)
#   pert <- seq(0.001, 0.09, length.out = 500)
#   logl.values <- c()
#   init.par    <- c() 
#   for (r in 1:R){
#     shift        <- seq(-pert[r], pert[r], length.out = Q)
#     lc.beta      <- lc.delta  <- lc.lnsigma_e  <- lc.lnsigma_v  <- lc.athrho  <- c()
#     for (i in 1:Q) {
#       lc.beta       <- c(lc.beta,      beta      + shift[i])
#       lc.delta      <- c(lc.delta,     delta     + shift[i])
#       lc.lnsigma_e  <- c(lc.lnsigma_e, lnsigma_e + shift[i])
#       lc.lnsigma_v  <- c(lc.lnsigma_v, lnsigma_v + shift[i])
#       lc.athrho     <- c(lc.athrho,    athrho    + shift[i])
#     }
#     thetas <- c(lc.beta, lc.delta, s_class, lc.lnsigma_e, lc.lnsigma_v, lc.athrho)
#     init.par <- rbind(init.par, thetas)
#     logl.values <- c(logl.values, 
#                      sum(ml_lcivd2(thetas, 
#                                    y1 = y1, 
#                                    y2 = y2, 
#                                    X  = X, 
#                                    Z  = Z, 
#                                    W  = Hl, 
#                                    Q = Q, 
#                                    gradient = FALSE)))
#   }
#   start <- init.par[which(max(logl.values) == logl.values), ]
#   names(beta)      <- paste("eq.1", colnames(X), sep = ".")
#   names(delta)     <- paste("eq.2", colnames(Z), sep = ".")
#   names(lnsigma_e) <- paste("eq.1", "lnsigma", sep = ".")
#   names(lnsigma_v) <- paste("eq.2", "lnsigma", sep = ".")
#   lc.nbeta         <- lc.ndelta  <- lc.nlnsigma_e  <- lc.nlnsigma_v  <- lc.nathrho  <- c()
#   for (i in 1:Q) {
#     lc.nbeta      <- c(lc.nbeta ,     paste('class', i, names(beta)      , sep = '.'))
#     lc.ndelta     <- c(lc.ndelta,     paste('class', i, names(delta)     , sep = '.'))
#     lc.nlnsigma_e <- c(lc.nlnsigma_e, paste('class', i, names(lnsigma_e) , sep = '.'))
#     lc.nlnsigma_v <- c(lc.nlnsigma_v, paste('class', i, names(lnsigma_v) , sep = '.'))
#     lc.nathrho    <- c(lc.nathrho,    paste('class', i, "athrho"         , sep = '.'))
#   }
#   names(lc.beta)      <- lc.nbeta
#   names(lc.delta)     <- lc.ndelta
#   names(lc.lnsigma_e) <- lc.nlnsigma_e
#   names(lc.lnsigma_v) <- lc.nlnsigma_v
#   names(lc.athrho)    <- lc.nathrho
#   names(start) <- c(lc.nbeta, lc.ndelta, namesH, lc.nlnsigma_e, lc.nlnsigma_v, lc.nathrho)
#   out <- list(start = start, pert.init = pert[which(max(logl.values) == logl.values)])
#   return(out)
# }


# ml_lcivd <- function(theta, y1, y2, X, Z, W, Q, gradient = TRUE){
#   #param: bq, dq, lq, lnsigma_eq, lnsigma_vq, athrho_q
#   K <- ncol(X)
#   P <- ncol(Z)
#   L <- ncol(W[[1]])
#   N <- nrow(X)
#   
#   beta      <- matrix(theta[1L:(K * Q)], nrow = K, ncol = Q)                        # Matrix of K*Q
#   delta     <- matrix(theta[((K * Q) + 1):((P * Q) + (K * Q))], nrow = P, ncol = Q) # Matrix of P*Q
#   lambda    <- theta[((P * Q) + (K * Q) + 1):(L + (P * Q) + (K * Q))]               # vector of (Q - 1) * L
#   lnsigma_e <- repRows(tail(theta, n = 3 * Q)[1:Q], n = nrow(X))                    # Matrix of n * Q
#   lnsigma_v <- repRows(tail(theta, n = 2 * Q)[1:Q], n = nrow(X))                    # Matrix of n * Q
#   athrho    <- repRows(tail(theta, n = Q), n = nrow(X))                             # Matrix of n * Q
#   sigma_e   <- exp(lnsigma_e)
#   sigma_v   <- exp(lnsigma_v)
#   rho       <- tanh(athrho)
#   
#   rownames(beta)  <- colnames(X); colnames(beta)  <- paste("class", 1:Q, sep = ":")
#   rownames(delta) <- colnames(Z); colnames(delta) <- paste("class", 1:Q, sep = ":")
#   
#   # Create components for log-likelihood function 
#   f      <- dnorm
#   index1 <- tcrossprod(X, t(beta))        #N * Q
#   index2 <- tcrossprod(Z, t(delta))       #N * Q
#   y1     <- repCols(y1, n = Q)            #N * Q
#   y2     <- repCols(y2, n = Q)            #N * Q
#   a      <- (y1 - index1)
#   b      <- (y2 - index2)
#   r      <- (sigma_e / sigma_v)
#   tiq    <- (a - r * rho * b)
#   aiq    <- - (1 / (2 * (1 - rho^2) * sigma_e^2)) * (tiq) ^ 2
#   biq    <- b / sigma_v
#   #P1     <- (1 / (sqrt(2 * pi * (1 - rho^2) * sigma_e^2))) * exp(aiq)
#   P1     <- (1 / (sqrt(2 * pi * (1 - rho^2) * sigma_e^2))) * exp(aiq)
#   #P1 <-  (1 / sqrt((1 - rho^2) * sigma_e^2)) * f(tiq / sqrt((1 - rho^2) * sigma_e^2))
#   P1     <- pmax(P1, .Machine$double.eps)
#   P1[is.na(P1)] <- 0
#   P2   <- (1 / sigma_v) * f(biq)
#   P1[is.na(P2)] <- 0
#   
#   # Make weights from a multinomial logit model
#   ew <- lapply(W, function(x) exp(crossprod(t(x), lambda)))
#   sew <- suml(ew)
#   Wiq <- lapply(ew, function(x){ v <- x / sew;
#   v[is.na(v)] <- 0;
#   as.vector(v)})
#   Wiq <- Reduce(cbind, Wiq) 
#   
#   Piq <- P1 * P2
#   #Piq <- pmax(P1 * P2 , .Machine$double.eps)
#   Pi <- pmax(apply(Wiq * Piq, 1, sum), .Machine$double.eps)
#   Pi[is.na(Pi)] <- 0
#   LL <- log(Pi)
#   
#   ##Gradient
#   if (gradient) {
#     ff       <- function(x) -x * dnorm(x)
#     #frat     <- function(x) ff(x) / f(x)
#     dp1  <- tiq / ((1 - rho^2) * sigma_e^2)
#     #dp1 <- -1 * (ff(tiq/sqrt((1 - rho^2) * sigma_e^2)) / f(tiq/sqrt((1 - rho^2) * sigma_e^2))) * (1/sqrt((1 - rho^2) * sigma_e^2))
#     dp2  <- -1 * (dp1 * r * rho + (ff(biq) / f(biq)) * (1 / sigma_v))
#     #dp2   <- (ff(tiq/sqrt((1 - rho^2) * sigma_e^2)) / f(tiq/sqrt((1 - rho^2) * sigma_e^2))) * (1/sqrt((1 - rho^2) * sigma_e^2)) * r * rho - (ff(biq) / f(biq)) * (1 / sigma_v)
#     dpse <- -1 + (a * tiq) / ((1 -  rho^2) * sigma_e^2)
#     dpsv <- - (tiq * rho * b) / ((1 -rho^2) * sigma_e * sigma_v)  - (ff(biq) / f(biq) * b / sigma_v) - 1
#     dpdt <-   (rho * sech(athrho)^2) / (1 - rho^2) - ((sech(athrho)^2 * tiq) / ((1 -rho^2)^2 * sigma_e^2)) * (a * rho - b*r)
#     
#     Qiq  <-  (Wiq * Piq / Pi)
#     eta  <-  Qiq * dp1  
#     etar <- eta[, rep(1:Q, each = K)]
#     Xg <- X[, rep(1:K, Q)] # N * (K * Q)
#     grad.beta <- Xg * etar
#     
#     eta <- Qiq * dp2
#     etar <- eta[, rep(1:Q, each = P)]
#     Zg <- Z[, rep(1:P, Q)] # N * (P * Q)
#     grad.delta <- Zg * etar
#     
#     grad.lnsigma_e <- Qiq * dpse
#     grad.lnsigma_v <- Qiq * dpsv
#     
#     grad.athrho <- Qiq * dpdt
#     
#     Wg <- vector(mode = "list", length = Q)
#     IQ <- diag(Q)
#     for (q in 1:Q) Wg[[q]] <- rowSums(Qiq * (repRows(IQ[q, ], N) - repCols(Wiq[, q], Q)))
#     grad.lambda <- suml(mapply("*", W, Wg, SIMPLIFY = FALSE))
#     grad <- cbind(grad.beta, grad.delta, grad.lambda, grad.lnsigma_e, grad.lnsigma_v, grad.athrho)
#     colnames(grad) <- names(theta)
#     attr(LL,'gradient') <- grad
#   }
#   LL
# }

