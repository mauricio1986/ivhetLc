
#' @title Estimation of instrumental variables models with latent classes
#' @author Mauricio Sarrias
#' @import stats methods Formula maxLik utils flexmix
#' @export
ivhetLc <- function(formula, data, subset, na.action,
                     method     = "bfgs",
                     type.init  = c("lc", "iv"),
                     Q          = 2,
                     print.init = FALSE,
                     messages   = TRUE,
                     start      =  NULL,
                     pert.init  = 0.05,
                     ...){

  callT     <- match.call(expand.dots = TRUE)
  callF     <- match.call(expand.dots = FALSE)
  nframe    <- length(sys.calls())
  type.init <- match.arg(type.init)

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
  if (length(f)[2L] < 2L) stop("You need to specify the instruments in the second part of the formula")


  # ============================
  # 4. Make variables
  # ============================

  # Extract variables for each part
  y1    <- model.response(mf)                                 # Continuous dependent variable
  y.var <- f[[2]]                                             # Name of continuous dependent variable
  X     <- model.matrix(f, data = mf, rhs = 1)                # Variables in the first equation
  Z     <- model.matrix(f, data = mf, rhs = 2)                # Variables in the second equation
  y2    <- X[, !(colnames(X) %in% colnames(Z)), drop = FALSE] # Endogenous continuous variable
  end.var     <- colnames(y2)                                 # Name of endogenous variable
  instruments <- colnames(Z)                                  # Name of instruments
  K <- ncol(X)
  P <- ncol(Z)
  # Checks
  if (P < K) stop("Model underidentified: check your variables in formula")
  if (messages && (K == P)) cat("\nEstimating a just identified model....\n")
  if (messages && (P > K))  cat("\nEstimating an overidentified model....\n")
  if (length(end.var) > 1L) stop(" ivhetLc only works with one endogenous variable")


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
    if (type.init == "iv"){
      s_class <- rep(0, ncol(H))
      names(s_class) <- namesH
      if(messages) cat("Estimating an IVs model for initial values", fill = TRUE)
      #ols1         <- lm(y1 ~ X - 1)
      iv           <- AER::ivreg(y1 ~ X - 1| Z - 1)
      ols2         <- lm(y2 ~ Z - 1)
      beta         <- coef(iv)
      delta        <- coef(ols2)
      lnsigma_e    <- log(iv$sigma)
      lnsigma_v    <- log(sigma(ols2))
      rho          <- cor(resid(iv), resid(ols2))
      athrho       <- 0.5 * log((1 + rho)/(1 - rho))
      names(beta)  <- paste("eq.1", colnames(X), sep = ".")
      names(delta) <- paste("eq.2", colnames(Z), sep = ".")
      names(lnsigma_e) <- paste("eq.1", "lnsigma", sep = ".")
      names(lnsigma_v) <- paste("eq.2", "lnsigma", sep = ".")
      init.shift   <- cumprod(c(1, rep(pert.init + 1, Q - 1)))
      lc.beta      <- lc.delta  <- lc.lnsigma_e  <- lc.lnsigma_v  <- lc.athrho  <- c()
      lc.nbeta     <- lc.ndelta  <- lc.nlnsigma_e  <- lc.nlnsigma_v  <- lc.nathrho  <- c()
      for (i in 1:Q) {
        lc.beta       <- c(lc.beta,      beta      * init.shift[i])
        lc.delta      <- c(lc.delta,     delta     * init.shift[i])
        lc.lnsigma_e  <- c(lc.lnsigma_e, lnsigma_e * init.shift[i])
        lc.lnsigma_v  <- c(lc.lnsigma_v, lnsigma_v * init.shift[i])
        lc.athrho     <- c(lc.athrho,    athrho    * init.shift[i])
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
    } else {
      if(messages) cat("Estimating a LC model for initial values", fill = TRUE)
      for1  <- formula(f, rhs = 1)
      end   <- colnames(X)[!(colnames(X) %in% colnames(Z))]
      for2  <- as.formula(paste(end, Reduce(paste, deparse(formula(f, lhs = 0, rhs = 2)))))
      for3  <- formula(f, lhs = 0, rhs = 3)
      lcm <- refit(flexmix(~., k = Q, data = mf,
                     model = list(FLXMRglm(for1),
                                  FLXMRglm(for2)),
                     concomitant = FLXPmultinom(for3)))
      s.par     <- lcm@coef[1:((ncol(X) + 1) *Q)]
      r.par     <- lcm@coef[(((ncol(X) + 1) *Q) + 1):((ncol(X) + ncol(Z) + 2)*Q)]
      lc.class  <- lcm@coef[c(-(1:((ncol(X) + ncol(Z) + 2) * Q)))]
      lc.athrho <- rep(0, Q)
      lc.beta   <- as.numeric(matrix(s.par, ncol = Q)[1:ncol(X), ])
      lc.delta  <- as.numeric(matrix(r.par, ncol = Q)[1:ncol(Z), ])
      lc.lnsigma_e <- as.numeric(matrix(s.par, ncol = Q)[(ncol(X) + 1), ])
      lc.lnsigma_v <- as.numeric(matrix(r.par, ncol = Q)[(ncol(Z) + 1), ])
      names.beta  <- paste("eq.1", colnames(X), sep = ".")
      names.delta <- paste("eq.2", colnames(Z), sep = ".")
      names.lnsigma_e <- paste("eq.1", "lnsigma", sep = ".")
      names.lnsigma_v <- paste("eq.2", "lnsigma", sep = ".")
      lc.nbeta     <- lc.ndelta  <- lc.nlnsigma_e  <- lc.nlnsigma_v  <- lc.nathrho  <- c()
      for (i in 1:Q) {
        lc.nbeta      <- c(lc.nbeta ,     paste('class', i, names.beta      , sep = '.'))
        lc.ndelta     <- c(lc.ndelta,     paste('class', i, names.delta     , sep = '.'))
        lc.nlnsigma_e <- c(lc.nlnsigma_e, paste('class', i, names.lnsigma_e , sep = '.'))
        lc.nlnsigma_v <- c(lc.nlnsigma_v, paste('class', i, names.lnsigma_v , sep = '.'))
        lc.nathrho    <- c(lc.nathrho,    paste('class', i, "athrho"         , sep = '.'))
      }
      start <- c(lc.beta, lc.delta, lc.class, lc.lnsigma_e, lc.lnsigma_v, lc.athrho)
      names(start) <- c(lc.nbeta, lc.ndelta, namesH, lc.nlnsigma_e, lc.nlnsigma_v, lc.nathrho)
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
  m <- match(c('method', 'print.level', 'iterlim',
               'start','tol', 'ftol', 'steptol', 'fixed', 'constraints',
               'control', 'gradient', 'hessian'),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]] <- as.name("maxLik")
  if(messages) cat("Estimating an IVLC-Linear model", "\n")
  opt$logLik <- as.name("ml_lcivd")
  opt[c("y1", "y2", "X", "Z", "W", "Q")] <- list(as.name("y1"),
                                                     as.name("y2"),
                                                     as.name("X"),
                                                     as.name("Z"),
                                                     as.name("Hl"),
                                                     as.name("Q"))
  out <- eval(opt, sys.frame(which = nframe))

  ##############################
  ## 6 Save results
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
  out$Q           <- Q
  out$y.var       <- y.var
  class(out) <- c("lciv", "maxLik", class(out))
  return(out)
}



ml_lcivd <- function(theta, y1, y2, X, Z, W, Q, gradient = TRUE, hessian = FALSE){
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
  #if (any(rho == 1 || rho == -1)) return(NA)
  if (any(rho == 1 | rho == -1)) return(NA)

  rownames(beta)  <- colnames(X); colnames(beta)  <- paste("class", 1:Q, sep = ":")
  rownames(delta) <- colnames(Z); colnames(delta) <- paste("class", 1:Q, sep = ":")

  # Create components for log-likelihood function
  index1 <- tcrossprod(X, t(beta))        #N * Q
  index2 <- tcrossprod(Z, t(delta))       #N * Q
  onesq  <- matrix(rep(1, Q), nrow = Q)
  y1     <- kronecker(t(onesq), y1)
  y2     <- kronecker(t(onesq), y2)        #N * Q
  r      <- (sigma_e / sigma_v)
  sd.y1  <- sqrt((1 - rho^2) * sigma_e^2)
  P1     <- dnorm(y1, mean = index1 + r * rho * (y2 - index2), sd = sd.y1)
  P2     <- dnorm(y2, mean = index2, sd = sigma_v)

  # Make weights from a multinomial logit model
  ew  <- lapply(W, function(x) exp(crossprod(t(x), lambda)))
  sew <- suml(ew)
  Wiq <- lapply(ew, function(x){ v <- x / sew;
                                v[is.na(v)] <- 0;
                                as.vector(v)})
  Wiq <- Reduce(cbind, Wiq)                          # N * Q

  Piq <- P1 * P2
  Pi  <- apply(Wiq * Piq, 1, sum)
  Pi  <- pmax(Pi, .Machine$double.eps)
  LL  <- log(Pi)


  ##Gradient
  if (gradient){
    w_iq   <- (Wiq * Piq / Pi)                                # N * Q
    a_iq   <- (y1 - index1 - r * rho * (y2 - index2)) / sd.y1 # N * Q
    b_iq   <- (y2 - index2) / sigma_v                         # N * Q
    onesk  <- matrix(rep(1, K), nrow = K)                     # K * 1
    onesp  <- matrix(rep(1, P), nrow = P)                     # P * 1
    Xg     <- kronecker(t(onesq), X)
    Zg     <- kronecker(t(onesq), Z)

    # da/dzeta and db/dzeta
    da_b   <- - (1 / sd.y1)
    da_d   <- (1 / sd.y1) * (r * rho)
    da_lne <- - (y1 - index1) / sd.y1
    da_lnv <- ((sigma_e * rho) / sd.y1) * b_iq
    da_t   <- ((y1 - index1) * rho - b_iq * sigma_e) / sd.y1
    db_d   <- - (1/ sigma_v)
    db_lnv <- - b_iq

    # dlnp / dzeta
    dlp_b   <- Xg * kronecker(-a_iq * da_b, t(onesk))
    dlp_d   <- Zg * kronecker(-a_iq * da_d - b_iq * db_d, t(onesp))
    dlp_lne <- -1 - a_iq * da_lne
    dlp_lnv <- -1 - a_iq * da_lnv - b_iq * db_lnv
    dlp_t   <- rho - a_iq * da_t

    # g.zeta
    g_b   <- dlp_b * kronecker(w_iq, t(onesk))
    g_d   <- dlp_d * kronecker(w_iq, t(onesp))
    g_lne <- w_iq * dlp_lne
    g_lnv <- w_iq * dlp_lnv
    g_t   <- w_iq * dlp_t

    # dlpi / dlambda
    Wg <- vector(mode = "list", length = Q)
    IQ <- diag(Q)
    for (q in 1:Q) Wg[[q]] <- rowSums(w_iq * (repRows(IQ[q, ], N) - repCols(Wiq[, q], Q))) # N * Q to N * 1
    g_lambda <- suml(mapply("*", W, Wg, SIMPLIFY = FALSE))


    g.theta <- cbind(g_b, g_d, g_lambda, g_lne, g_lnv, g_t)
    colnames(g.theta) <- names(theta)
    g.theta[is.nan(g.theta)] <- 0
    attr(LL,'gradient') <- g.theta

  }
  if (hessian){
    # Construct own hessian
    H <- matrix(0, (K + P + 3) * Q + L * (Q - 1), (K + P + 3) * Q + L * (Q - 1))
    colnames(H) <- rownames(H) <- names(theta)
    D  <- crossprod(g.theta, g.theta)
    dlnp <- cbind(dlp_b,
                  dlp_d,
                  suml(mapply("*", W, Wg, SIMPLIFY = FALSE)),
                  dlp_lne,
                  dlp_lnv,
                  dlp_t)
    wD <- crossprod(g.theta, dlnp)
    for (q in 1:Q){
      wH_bb     <- - crossprod(w_iq[, q] * X * da_b[, q], X * da_b[, q]) # K * K
      wH_bd     <- - crossprod(w_iq[, q] * X * da_b[, q], Z * da_d[, q]) # K * P
      wH_blne   <- colSums(- w_iq[, q] * (X * da_b[, q] * da_lne[, q] + a_iq[, q] * (- X * da_b[, q]))) # 1 * K
      wH_blnv   <- colSums(-w_iq[, q] * X * da_b[, q] * da_lnv[, q])
      wH_bt     <- colSums(-w_iq[, q]*(X * da_b[, q] * da_t[, q] + a_iq[, q] * (- X * (rho[, q] / sigma_v[, q])*(1/sd.y1[, q]))))
      wH_dd     <- - crossprod(w_iq[, q] * Z * da_d[, q], Z * da_d[, q]) - crossprod(w_iq[, q] * Z * db_d[, q], Z * db_d[, q])
      wH_dlne   <- colSums(- w_iq[, q] * Z * da_d[, q] * da_lne[, q])
      wH_dlnv   <- colSums(- w_iq[, q] * (Z * da_d[, q] * da_lnv[, q] + a_iq[, q] * (-Z * da_d[, q]))) +
                   colSums(- w_iq[, q] * (Z * db_d[, q] * db_lnv[, q] + b_iq[, q] * (-Z * db_d[, q])))
      wH_dt     <- colSums(- w_iq[, q] * (Z * da_d[, q] * da_t[, q] + a_iq[, q] * (Z * (1 / sigma_v[, q]) * (1 / sqrt(1 - rho[, q]^2)))))
      wH_lnelne <- sum(-w_iq[, q] * (da_lne[, q]^2 + a_iq[, q]* (-da_lne[, q])))
      wH_lnelnv <- sum(-w_iq[, q] * da_lne[, q] * da_lnv[, q])
      wH_lnet   <- sum(-w_iq[, q] * (da_lne[, q] * da_t[, q]  + a_iq[, q] * (da_lne[, q] * rho[, q])))
      wH_lnvlnv <- sum(-w_iq[, q] * (da_lnv[, q]^2 + a_iq[, q] * (-da_lnv[, q]))) +
                   sum(-w_iq[, q] * (db_lnv[, q]^2 + b_iq[, q] * (-db_lnv[, q])))
      wH_lnvt   <- sum(-w_iq[, q] * (da_lnv[, q] * da_lnv[, q] + a_iq[, q] * (-b_iq[, q]/sd.y1[, q])))
      wH_tt     <- sum(-w_iq[, q] * ((1 -  rho[, q]) - (da_t[, q]^2 + a_iq[, q] * (((y1[, q] - index1[, q]) - b_iq[, q] * sigma_e[, q] * rho[, q]) / sd.y1[, q]))))
      ik        <- c((q - 1) * K + 1):(q * K)
      ip        <- c(K*Q + (q - 1) * P + 1):(K*Q + q*P)
      ic        <- (K + P) * Q + L * (Q - 1) + q
      H[ik, ik]       <- wH_bb + D[ik, ik] + wD[ik, ik]
      H[ik, ip]       <- wH_bd + D[ik, ip] + wD[ik, ip]
      H[ik, ic]       <- wH_blne + D[ik, ic] + wD[ik, ic]
      H[ik, ic + Q]     <- wH_blnv + D[ik, ic + Q]   + wD[ik, ic + Q]
      H[ik, ic + 2*Q]   <- wH_bt   + D[ik, ic + 2*Q] + wD[ik, ic + 2*Q]
      H[ip, ik]         <- t(H[ik, ip])
      H[ip, ip]         <- wH_dd   + D[ip, ip] + wD[ip, ip]
      H[ip, ic]         <- wH_dlne + D[ip, ic] + wD[ip, ic]
      H[ip, ic + Q]     <- wH_dlnv + D[ip, ic + Q] + wD[ip, ic + Q]
      H[ip, ic + 2*Q]   <- wH_dt   + D[ip, ic + 2*Q] + wD[ip, ic + 2*Q]
      H[ic, ik]         <- H[ik, ic]
      H[ic, ip]         <- H[ip, ic]
      H[ic, ic]         <- wH_lnelne + D[ic, ic] + wD[ic, ic]
      H[ic, ic + Q]     <- wH_lnelnv + D[ic, ic + Q] + wD[ic, ic + Q]
      H[ic, ic + 2*Q]   <- wH_lnet + D[ic, ic + 2*Q] + wD[ic, ic + 2*Q]
      H[ic + Q, ik]     <- t(H[ik, ic + Q])
      H[ic + Q, ip]     <- t(H[ip, ic + Q])
      H[ic + Q, ic]     <- H[ic, ic + Q]
      H[ic + Q, ic + Q] <- wH_lnvlnv + D[ic + Q, ic + Q] + wD[ic + Q, ic + Q]
      H[ic + Q, ic + 2*Q] <- wH_lnvt + D[ic + Q, ic + 2 * Q] + wD[ic + Q, ic + 2 * Q]
      H[ic + 2*Q, ik]         <- t(H[ik, ic + 2*Q])
      H[ic + 2*Q, ip]         <- t(H[ip, ic + 2*Q])
      H[ic + 2*Q, ic] <- H[ic, ic + 2*Q]
      H[ic + 2 * Q, ic + Q] <- H[ic + Q, ic + 2*Q]
      H[ic + 2* Q, ic + 2*Q] <- wH_tt + D[ic + 2 * Q, ic + 2 * Q] + wD[ic + 2 * Q, ic + 2 * Q]
    }
    #dlp_dbb  <- - da_db^2
    #dlp_dbd  <- - da_db * da_dd
    #dlp_dblne <-  -(da_db * da_lne + a_iq * (-da_db))
    #dlp_dblnv <- - (da_db * da_lnv)
    #dlp_dbt   <- - (da_db * da_t + a_iq * )
    #dp_zeta_q   <- cbind(grad.beta, grad.delta, grad.lnsigma_e, grad.lnsigma_v, grad.athrho)
    #w_dp_zeta_q <- t(dp_zeta_q) %*% dp_zeta_q
    print(H)
  }
  attr(LL, 'pi')      <- Wiq
  attr(LL, 'post_pi') <- Wiq * Piq / Pi
  LL
}





